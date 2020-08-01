from cc3d.core.PySteppables import *

import math
import random

# Conversion Factors
s_to_mcs = 5 * 60  # s/mcs

# =============================
# CompuCell Parameters
# cell
cell_diameter = 5.0
cell_volume = cell_diameter ** 2
volume_lm = 9

# virus diffusion
virus_dc = 1.875

# cytokine diffusion
cytokine_dc = 7.5

# Lambda Chemotaxis
lamda_chemotaxis = 1E5

# Antimony/SBML model step size
sbml_step_size = s_to_mcs / 60 / 60 / 24

# Initial fraction of infected cells
init_infect = 0.05

# Option to use local or ODE virus
use_local_virus = True

ode_model_string = f'''model ODEModel()
    //Equations
    E1: T -> I1 ; beta*T*V ; //Infection rate
    E2: I1 -> I2 ; k*I1 ; //Infection rate
    E3: -> V ; p*I2 ; //Virus Production
    E4: V -> ; c*V ; // Virus Decay
    E5: I2 -> D ; d*I2 ; // Infected Cell Clearance (apopotosis + innate immune response)
    E6: -> C ; pc*(I1+I2) ; // Cytokine production
    E7: C -> ; cc*C ; // Cytokine decay
    E8: C -> Cl ; kc * C // Cytokine transport to the lymph node
    E9: Cl -> ; ccl * Cl // Lymph node cytokine decay
    E10: -> El ; pel*Cl + rel*Kel*Cl*El/(Kel + El) ;
    E11: El -> E; ke*El   ;
    E12: E ->  ; dE*E     ;
    E13: I2 -> D ; kei2*dei2*E/(E+kei2+I2)*I2 ;


    //Parameters
    beta = 6.2E-5 ; // 1.0/(TCID*day)
    k = 4.0 ; // 1.0/day
    p = 1.0 ; // TCID/cell/day
    c = 9.4 ; // 1.0/day
    d = 2.4E-1; // 1.0/day
    pc = 1.0 ;
    cc = 2.0 ; // 1.0/day
    kc = 0.5 ; // 1.0/day
    ccl = 0.5 ; // 1.0/day
    pel = 1E-4 ;
    rel = 0.005 ;
    Kel = 100 ;
    ke = 0.1 ;
    dE = 0.5 ;
    dei2 = 1E-2 ; //15E3
    kei2 = 5E2 ;

    //Killing mechanisms
    -> viralDeath ; d*I2
    -> cd8Death ; kei2*dei2*E/(E+kei2+I2)*I2
    viralDeath = 0 ;
    cd8Death = 0 ; 

    //Inputs
    V = 0.0

    //Initial Conditions
    T0 = 1E7
    T = T0
    I1 = 75 
    C = 0.0
    E = 0.0
    El = 0.0
    end'''

lymph_node_model_string = f'''model LymphNodeMode()
    //Equations
    E8: -> Cl ; kc * C // Cytokine transport to the lymph node
    E9: Cl -> ; ccl * Cl // Lymph node cytokine decay
    E10: -> El ; pel*Cl + rel*Kel*Cl*El/(Kel + El) ;
    E11: El -> ; ke*El   ;

    //Parameters
    beta = 6.2E-5 ; // 1.0/(TCID*day)
    k = 4.0 ; // 1.0/day
    p = 1.0 ; // TCID/cell/day
    c = 9.4 ; // 1.0/day
    d = 2.4E-1; // 1.0/day
    pc = 1.0 ;
    cc = 2.0 ; // 1.0/day
    kc = 0.5 ; // 1.0/day
    ccl = 0.5 ; // 1.0/day
    pel = 1E-4 ;
    rel = 0.005 ;
    Kel = 100 ;
    ke = 0.1 ;
    dE = 0.5 ;
    dei2 = 1E-2 ; //15E3
    kei2 = 5E2 ;
    T0 = 1E7 ;

    //Inputs
    C = 0.0

    //Outputs
    El = 0.0
    end'''


class ModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to ODE model solver
        self.rr = None

        # Reference to ODE model solver without spatialization
        self.rr_ode = None

        # Population data window and path
        self.pop_data_win = None

        # Diffusion data window and path
        self.med_diff_data_win = None

        # Death data window and path
        self.death_data_win = None

        # Cell death mechanism tracking
        self.death_mech = {'viral': 0,
                           'contact': 0}

        # Scaling coefficients
        self.scale_by_volume = None
        self.scale_by_voxel = None
        self.scale_by_time = None

    def start(self):

        # Step: Initialize epithelial sheet

        # Enforce compatible lattice dimensions with epithelial cell size
        assert self.dim.x % cell_diameter == 0 and self.dim.y % cell_diameter == 0, \
            f'Lattice dimensions must be multiples of the unitless cell diameter (currently cell_diameter = {cell_diameter})'

        for x in range(0, self.dim.x, int(cell_diameter)):
            for y in range(0, self.dim.y, int(cell_diameter)):
                cell = self.new_cell(self.UNINFECTED)
                self.cell_field[x:x + int(cell_diameter), y:y + int(cell_diameter), 0] = cell

                cell.targetVolume = cell_volume
                cell.lambdaVolume = volume_lm

        # Infect a cell
        num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))
        if init_infect == 0:
            cell = self.cell_field[self.dim.x // 2, self.dim.y // 2, 0]
            cell.type = self.INFECTED
        else:
            num_to_infect = int(init_infect * num_ecs)
            num_infected = 0
            while num_infected < num_to_infect:
                xi = random.randint(0, self.dim.x)
                yi = random.randint(0, self.dim.y)
                cell = self.cell_field[xi, yi, 0]
                if cell is not None and cell.type == self.UNINFECTED:
                    cell.type = self.INFECTED
                    num_infected += 1

        # Step: initialize ODEs

        #   Generate solver instance
        self.add_free_floating_antimony(model_string=lymph_node_model_string,
                                        model_name='LymphNodeMode',
                                        step_size=sbml_step_size)

        #   Get reference to solver
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == 'LymphNodeMode':
                self.rr = rr

        # Step: calculate scaling coefficients
        self.scale_by_volume = num_ecs / self.rr["T0"]
        self.scale_by_voxel = 1 / self.rr["T0"] / cell_volume
        self.scale_by_time = s_to_mcs / 60 / 60 / 24

        # Step: load secretion values from ODEs into XML

        #   Extracellular virus: diffusion and decay
        self.get_xml_element('virus_dc').cdata = virus_dc
        self.get_xml_element('virus_decay').cdata = self.rr["c"] * self.scale_by_time

        #   Local cytokine: diffusion and decay
        self.get_xml_element('cytokine_dc').cdata = cytokine_dc
        self.get_xml_element('cytokine_decay').cdata = (self.rr["cc"] + self.rr["kc"]) * self.scale_by_time

        # Step: data tracking setup

        # Initialize windows and paths if requested

        dot_size = 5
        line_size = 3

        # Population data tracking plot setup

        self.pop_data_win = self.add_new_plot_window(title='Population data',
                                                     x_axis_title='MCS',
                                                     y_axis_title='Numer of cells',
                                                     x_scale_type='linear',
                                                     y_scale_type='log',
                                                     grid=True,
                                                     config_options={'legend': True})

        self.pop_data_win.add_plot("Uninfected", style='Dots', color='blue', size=dot_size)
        self.pop_data_win.add_plot("Infected", style='Dots', color='red', size=dot_size)
        self.pop_data_win.add_plot("VirusReleasing", style='Dots', color='green', size=dot_size)
        self.pop_data_win.add_plot("Dead", style='Dots', color='yellow', size=dot_size)
        self.pop_data_win.add_plot("CD8Local", style='Dots', color='white', size=dot_size)
        self.pop_data_win.add_plot("CD8Lymph", style='Dots', color='purple', size=dot_size)

        self.pop_data_win.add_plot("UninfectedODE", style='Lines', color='blue', size=line_size)
        self.pop_data_win.add_plot("InfectedODE", style='Lines', color='red', size=line_size)
        self.pop_data_win.add_plot("VirusReleasingODE", style='Lines', color='green', size=line_size)
        self.pop_data_win.add_plot("DeadODE", style='Lines', color='yellow', size=line_size)
        self.pop_data_win.add_plot("CD8LocalODE", style='Lines', color='white', size=line_size)
        self.pop_data_win.add_plot("CD8LymphODE", style='Lines', color='purple', size=line_size)

        # Diffusive field data tracking setup

        self.med_diff_data_win = self.add_new_plot_window(title='Total diffusive species',
                                                          x_axis_title='MCS',
                                                          y_axis_title='Number of diffusive species per volume',
                                                          x_scale_type='linear',
                                                          y_scale_type='log',
                                                          grid=True,
                                                          config_options={'legend': True})

        self.med_diff_data_win.add_plot("MedViral", style='Dots', color='red', size=dot_size)
        self.med_diff_data_win.add_plot("MedCytLocal", style='Dots', color='blue', size=dot_size)
        self.med_diff_data_win.add_plot("MedCytLymph", style='Dots', color='green', size=dot_size)

        self.med_diff_data_win.add_plot("MedViralODE", style='Lines', color='red', size=line_size)
        self.med_diff_data_win.add_plot("MedCytLocalODE", style='Lines', color='blue', size=line_size)
        self.med_diff_data_win.add_plot("MedCytLymphODE", style='Lines', color='green', size=line_size)

        # Death mechanism data tracking setup

        self.death_data_win = self.add_new_plot_window(title='Death data',
                                                       x_axis_title='MCS',
                                                       y_axis_title='Numer of cells',
                                                       x_scale_type='linear',
                                                       y_scale_type='log',
                                                       grid=True,
                                                       config_options={'legend': True})

        self.death_data_win.add_plot("Viral", style='Dots', color='blue', size=5)
        self.death_data_win.add_plot("Contact", style='Dots', color='green', size=5)

        self.death_data_win.add_plot("ViralODE", style='Lines', color='blue', size=line_size)
        self.death_data_win.add_plot("ContactODE", style='Lines', color='green', size=line_size)

        #   Generate solver instance
        self.add_free_floating_antimony(model_string=ode_model_string,
                                        model_name='ODEModel',
                                        step_size=sbml_step_size)

        #   Get reference to solver
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == 'ODEModel':
                self.rr_ode = rr

        #   Apply initial conditions
        self.rr_ode["T"] = len(self.cell_list_by_type(self.UNINFECTED)) / self.scale_by_volume
        self.rr_ode["I1"] = len(self.cell_list_by_type(self.INFECTED)) / self.scale_by_volume
        self.rr_ode["I2"] = len(self.cell_list_by_type(self.VIRUSRELEASING)) / self.scale_by_volume
        self.rr_ode["D"] = len(self.cell_list_by_type(self.DEAD)) / self.scale_by_volume

    def step(self, mcs):
        # Prep stuff
        num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))
        virus_secretor = self.get_field_secretor("Virus")
        cytokine_secretor = self.get_field_secretor("cytokine")

        # Step: Immune cell killing
        dei2 = self.rr["dei2"] * self.scale_by_time / self.scale_by_volume

        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            cd8_area = 0
            srf_area = 0
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor is None or neighbor.type == self.CD8LOCAL:
                    srf_area += common_surface_area
                if neighbor is not None and neighbor.type == self.CD8LOCAL:
                    cd8_area += common_surface_area

            if cd8_area == 0:
                continue

            death_rate = dei2 * num_ecs / srf_area * cd8_area
            pr_death = 1 - math.exp(-death_rate)
            if random.random() < pr_death:
                cell.type = self.DEAD
                cell.targetVolume = 0
                self.death_mech['contact'] += 1

        # Plot contact death data
        if self.death_mech['contact'] > 0:
            self.death_data_win.add_data_point("Contact", mcs, self.death_mech['contact'])

        # Step: Viral killing

        death_rate = self.rr["d"] * self.scale_by_time
        pr_death = 1 - math.exp(-death_rate)
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            if random.random() < pr_death:
                cell.type = self.DEAD
                self.death_mech['viral'] += 1

        # Plot viral death data
        if self.death_mech['viral'] > 0:
            self.death_data_win.add_data_point("Viral", mcs, self.death_mech['viral'])

        # Step: Viral life cycle

        # Infected -> Virus releasing
        release_rate = self.rr["k"] * self.scale_by_time
        pr_release = 1 - math.exp(-release_rate)
        for cell in self.cell_list_by_type(self.INFECTED):
            if random.random() < pr_release:
                cell.type = self.VIRUSRELEASING

        # Internalization
        if use_local_virus:
            beta = self.rr["beta"] * self.scale_by_time / self.scale_by_voxel  # * num_ecs
            for cell in self.cell_list_by_type(self.UNINFECTED):
                seen_amount = virus_secretor.amountSeenByCell(cell) / cell.volume * self.dim.z
                infect_rate = beta * seen_amount
                pr_infect = 1 - math.exp(-infect_rate)
                if random.random() < pr_infect:
                    cell.type = self.INFECTED
        else:
            beta = self.rr["beta"] * self.scale_by_time
            infect_rate = beta * self.rr_ode["V"]
            pr_infect = 1 - math.exp(-infect_rate)
            for cell in self.cell_list_by_type(self.UNINFECTED):
                if random.random() < pr_infect:
                    cell.type = self.INFECTED

        # Plot epithelial population data

        num_cells_uninfected = len(self.cell_list_by_type(self.UNINFECTED))
        num_cells_infected = len(self.cell_list_by_type(self.INFECTED))
        num_cells_virusreleasing = len(self.cell_list_by_type(self.VIRUSRELEASING))
        num_cells_dead = len(self.cell_list_by_type(self.DEAD))
        if num_cells_uninfected > 0:
            self.pop_data_win.add_data_point('Uninfected', mcs, num_cells_uninfected)
        if num_cells_infected > 0:
            self.pop_data_win.add_data_point('Infected', mcs, num_cells_infected)
        if num_cells_virusreleasing > 0:
            self.pop_data_win.add_data_point('VirusReleasing', mcs, num_cells_virusreleasing)
        if num_cells_dead > 0:
            self.pop_data_win.add_data_point('Dead', mcs, num_cells_dead)

        # Step: Immune recruitment

        # Pass spatial information to ODEs
        self.rr["C"] = cytokine_secretor.totalFieldIntegral() / self.dim.z / self.scale_by_volume

        # Integrate ODEs
        self.rr.timestep()

        # Calculate outflow
        remove_rate = self.rr["dE"] * self.scale_by_time
        pr_remove = 1 - math.exp(-remove_rate)
        for cell in self.cell_list_by_type(self.CD8LOCAL):
            if random.random() < pr_remove:
                cell.targetVolume = 0

        # Calculate inflow
        num_add = 0
        add_rate = self.rr["El"] * self.rr["ke"] * self.scale_by_time * self.scale_by_volume
        exp_term = math.exp(-add_rate)
        sum_term = 1.0
        while random.random() < 1 - exp_term * sum_term:
            num_add += 1.0
            sum_term += add_rate ** num_add / math.factorial(num_add)

        for _ in range(int(num_add)):
            sample_frac = 0.01  # Number of sites to sample
            n_sites = self.dim.x * self.dim.y
            n_sites_frac = int(n_sites * sample_frac)

            max_concentration = 0
            x_seed = None
            y_seed = None
            sites_sampled = 0

            med_pixel_set = [ptd.pixel for ptd in self.pixel_tracker_plugin.getMediumPixelSet()]
            random.shuffle(med_pixel_set)

            for pixel in med_pixel_set:
                xi = pixel.x
                yi = pixel.y

                # No placing immune cells over a boundary
                if xi + int(cell_diameter / 2) >= self.dim.x or yi + int(cell_diameter / 2) >= self.dim.y:
                    continue

                open_space = True
                for x in range(xi, xi + int(cell_diameter / 2)):
                    for y in range(yi, yi + int(cell_diameter / 2)):
                        if self.cell_field[x, y, 1]:
                            open_space = False
                            break
                if open_space:
                    concentration_iteration = self.field.cytokine[xi, yi, 1]
                    sites_sampled += 1
                    if concentration_iteration >= max_concentration:
                        max_concentration = concentration_iteration
                        x_seed = xi
                        y_seed = yi
                if sites_sampled == n_sites_frac:
                    break
            if x_seed is not None:
                cell = self.new_cell(self.CD8LOCAL)
                self.cell_field[x_seed:x_seed + int(cell_diameter / 2), y_seed:y_seed + int(cell_diameter / 2), 1] = \
                    cell
                cell.targetVolume = cell_volume
                cell.lambdaVolume = volume_lm

        # Plot immune cell population data
        num_cells_immune = len(self.cell_list_by_type(self.CD8LOCAL))
        num_cells_immune_l = self.rr_ode["El"] * self.scale_by_volume
        if num_cells_immune > 0:
            self.pop_data_win.add_data_point('CD8Local', mcs, num_cells_immune)
        if num_cells_immune_l > 0:
            self.pop_data_win.add_data_point('CD8Lymph', mcs, num_cells_immune_l)

        # Step: Update chemotaxis

        for cell in self.cell_list_by_type(self.CD8LOCAL):
            cd = self.chemotaxisPlugin.getChemotaxisData(cell, "cytokine")
            if cd is None:
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "cytokine")
                cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
            concentration = self.field.cytokine[cell.xCOM, cell.yCOM, 1]
            cd.setLambda(lamda_chemotaxis / (1.0 + concentration))

        # Step: Secretion

        #   Extracellular virus: secretion by virus-releasing cells
        secr_amount = self.rr["p"] * self.scale_by_time
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            virus_secretor.secreteInsideCellTotalCount(cell, secr_amount / cell.volume)

        #   Local cytokine: secretion by infected and virus-releasing cells
        secr_amount = self.rr["pc"] * self.scale_by_time
        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING):
            cytokine_secretor.secreteInsideCellTotalCount(cell, secr_amount / cell.volume)

        # Plot total diffusive viral amount
        med_viral_total = virus_secretor.totalFieldIntegral() / self.dim.z
        med_cyt_total = cytokine_secretor.totalFieldIntegral() / self.dim.z
        med_cyt_lymph = self.rr["Cl"] * self.scale_by_volume
        if med_viral_total > 0:
            self.med_diff_data_win.add_data_point("MedViral", mcs, med_viral_total)
        if med_cyt_total > 0:
            self.med_diff_data_win.add_data_point("MedCytLocal", mcs, med_cyt_total)
        if med_cyt_lymph > 0:
            self.med_diff_data_win.add_data_point("MedCytLymph", mcs, med_cyt_lymph)

        # Step: Update ODE data and tracking

        self.rr_ode.timestep()

        # Population data tracking

        num_cells_uninfected = self.rr_ode["T"] * self.scale_by_volume
        num_cells_infected = self.rr_ode["I1"] * self.scale_by_volume
        num_cells_virusreleasing = self.rr_ode["I2"] * self.scale_by_volume
        num_cells_dead = self.rr_ode["D"] * self.scale_by_volume
        num_cells_immune = self.rr_ode["E"] * self.scale_by_volume
        num_cells_immune_l = self.rr_ode["El"] * self.scale_by_volume

        min_thresh = 0.1
        if num_cells_uninfected > min_thresh:
            self.pop_data_win.add_data_point('UninfectedODE', mcs, num_cells_uninfected)
        if num_cells_infected > min_thresh:
            self.pop_data_win.add_data_point('InfectedODE', mcs, num_cells_infected)
        if num_cells_virusreleasing > min_thresh:
            self.pop_data_win.add_data_point('VirusReleasingODE', mcs, num_cells_virusreleasing)
        if num_cells_dead > min_thresh:
            self.pop_data_win.add_data_point('DeadODE', mcs, num_cells_dead)
        if num_cells_immune > min_thresh:
            self.pop_data_win.add_data_point('CD8LocalODE', mcs, num_cells_immune)
        if num_cells_immune_l > min_thresh:
            self.pop_data_win.add_data_point('CD8LymphODE', mcs, num_cells_immune_l)

        # Diffusive field data tracking

        med_viral_total = self.rr_ode["V"] * self.scale_by_volume
        med_cyt_total = self.rr_ode["C"] * self.scale_by_volume
        med_cyt_lymph = self.rr_ode["Cl"] * self.scale_by_volume

        if med_viral_total > 0:
            self.med_diff_data_win.add_data_point("MedViralODE", mcs, med_viral_total)
        if med_cyt_total > 0:
            self.med_diff_data_win.add_data_point("MedCytLocalODE", mcs, med_cyt_total)
        if med_cyt_lymph > 0:
            self.med_diff_data_win.add_data_point("MedCytLymphODE", mcs, med_cyt_lymph)

        # Death mechanism data tracking

        num_viral = self.rr_ode['viralDeath'] * self.scale_by_volume
        num_contact = self.rr_ode['cd8Death'] * self.scale_by_volume

        min_thresh = 0.1
        if num_viral > min_thresh:
            self.death_data_win.add_data_point("ViralODE", mcs, num_viral)
        if num_contact > min_thresh:
            self.death_data_win.add_data_point("ContactODE", mcs, num_contact)

