from cc3d.core.PySteppables import *

import math
import random

# Data control options
plot_pop_data_freq = 10  # Plot population data frequency (disable with 0)
write_pop_data_freq = 0  # Write population data to simulation directory frequency (disable with 0)
plot_med_diff_data_freq = 10  # Plot total diffusive field amount frequency (disable with 0)
write_med_diff_data_freq = 0  # Write total diffusive field amount frequency (disable with 0)
plot_spat_data_freq = 0  # Plot spatial data frequency (disable with 0)
write_spat_data_freq = 0  # Write spatial data to simulation directory frequency (disable with 0)
plot_death_data_freq = 10  # Plot death data frequency (disable with 0)
write_death_data_freq = 0  # Write death data to simulation directory frequency (disable with 0)

plot_ode_sol = True  # Set to True to instantiate ODE model and plot, for comparison with spatial model

# Conversion Factors
s_to_mcs = 5 * 60  # s/mcs
um_to_lat_width = 2.0  # um/lattice_length

# Experimental Parameters
exp_cell_diameter = 10.0  # um

exp_virus_dc = 0.025  # um^2/s

exp_cytokine_dc_cyto = 0.1  # um^2/s; estimated diffusion constant in cytoplasm (B)
# ^ the  /10 is not experimental; added because of relative small area and because virus D is (or was) slowed down

# =============================
# CompuCell Parameters
# cell
cell_diameter = exp_cell_diameter / um_to_lat_width
cell_volume = cell_diameter ** 2
volume_lm = 9

# virus diffusion
virus_dc = exp_virus_dc * s_to_mcs / (um_to_lat_width ** 2)  # virus diffusion constant

cytokine_dc = exp_cytokine_dc_cyto * s_to_mcs / (um_to_lat_width ** 2)  # CK diff cst

# Lambda Chemotaxis
lamda_chemotaxis = 1E5

# Antimony/SBML model step size
sbml_step_size = 1.0

# Initial fraction of infected cells
init_infect = 0.05


def immune_recruitment_model_string_original(time_conv=1.0):
    """
    String generator for original ODE model
    :param time_conv: conversion factor from ODE model time scale to simulation time scale
    :return: model string
    """
    model_string = f'''model TcellModelOriginal()
    s_t = {time_conv}

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
    beta = 6.2E-5 * s_t ; // 1.0/(TCID*day)
    k = 4.0 * s_t ; // 1.0/day
    p = 1.0 * s_t ; // TCID/cell/day
    c = 9.4 * s_t ; // 1.0/day
    d = 2.4E-1 * s_t; // 1.0/day
    pc = 1.0 * s_t ;
    cc = 2.0 * s_t ; // 1.0/day
    kc = 0.5 * s_t ; // 1.0/day
    ccl = 0.5 * s_t ; // 1.0/day
    pel = 1E-4 * s_t ;
    rel = 0.005 * s_t ;
    Kel = 100 ;
    ke = 0.1 * s_t ;
    dE = 0.5 * s_t ;
    dei2 = 1E-2 * s_t ; //15E3
    kei2 = 5E2 ;
    
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
    return model_string


def immune_recruitment_model_string_ode(num_ec=1E7, num_infect=75.0, time_conv=1.0):
    """
    String generator for ODE model scaled to simulation domain
    :param num_ec: number of epithelial cells
    :param num_infect: number of initially infected cells
    :param time_conv: conversion factor from ODE model time scale to simulation time scale
    :return: model string
    """
    model_string = f'''model TcellModelODE()
    // Scaling coefficients
    s_v = {num_ec} / 1E7;
    s_t = {time_conv};
    
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
    beta = 6.2E-5 * s_t / s_v ; // 1.0/(TCID*day)
    k = 4.0 * s_t ; // 1.0/day
    p = 1.0 * s_t ; // TCID/cell/day
    c = 9.4 * s_t ; // 1.0/day
    d = 2.4E-1 * s_t; // 1.0/day
    pc = 1.0 * s_t ;
    cc = 2.0 * s_t ; // 1.0/day
    kc = 0.5 * s_t ; // 1.0/day
    ccl = 0.5 * s_t ; // 1.0/day
    pel = 1E-4 * s_t ;
    rel = 0.005 * s_t / s_v ;
    Kel = 100 * s_v ;
    ke = 0.1 * s_t ;
    dE = 0.5 * s_t ;
    dei2 = 1E-2 * s_t / s_v ; //15E3
    kei2 = 5E2 * s_v ;
    
    //Inputs
    V = 0.0
    
    //Initial Conditions
    T0 = {num_ec}
    T = T0
    I1 = {num_infect}
    C = 0.0
    E = 0.0
    El = 0.0

    // Death mechanism tracking
    -> viralDeath ; d*I2 ;
    -> cd8Death ; kei2*dei2*E/(E+kei2+I2)*I2 ;
    viralDeath = 0;
    cd8Death = 0;

    end'''
    return model_string


def immune_recruitment_model_string(num_ec=1E7, sites_per_cell=1.0, time_conv=1.0):
    """
    String generator for cellularized ODE model
    :param num_ec: number of epithelial cells
    :param sites_per_cell: number of sites per epithelial cell
    :param time_conv: conversion factor from ODE model time scale to simulation time scale
    :return: model string
    """
    model_string = f'''model TcellModel()
    // Scaling coefficients
    s_v = {num_ec} / 1E7;
    s_l = 1 / 1E7 / {sites_per_cell};
    s_t = {time_conv};

    //Equations
    E8: -> Cl ; kc * C // Cytokine transport to the lymph node
    E9: Cl -> ; ccl * Cl // Lymph node cytokine decay
    E10: -> El ; pel*Cl + rel*Kel*Cl*El/(Kel + El) ;
    E11: El -> ; ke*El   ;
    
    //Parameters
    beta = 6.2E-5 * s_t / s_l ; // 1.0/(TCID*day)
    k = 4.0 * s_t ; // 1.0/day
    p = 1.0 * s_t ; // TCID/cell/day
    c = 9.4 * s_t ; // 1.0/day
    d = 2.4E-1 * s_t; // 1.0/day
    pc = 1.0 * s_t ;
    cc = 2.0 * s_t ; // 1.0/day
    kc = 0.5 * s_t ; // 1.0/day
    ccl = 0.5 * s_t ; // 1.0/day
    pel = 1E-4 * s_t ;
    rel = 0.005 * s_t / s_v ;
    Kel = 100 * s_v ;
    ke = 0.1 * s_t ;
    dE = 0.5 * s_t ;
    dei2 = 1E-2 * s_t / s_v ; //15E3
    kei2 = 5E2 * s_v ;
    
    //Inputs
    V = 0.0
    
    //Initial Conditions
    T0 = {num_ec}
    T = T0
    C = 0.0
    E = 0.0
    El = 0.0
    end'''
    return model_string


class ModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to ODE model solver
        self.rr = None

        # Reference to ODE model solver without spatialization
        self.rr_ode = None

        # Population data window and path
        self.pop_data_win = None
        self.pop_data_path = None

        # Diffusion data window and path
        self.med_diff_data_win = None
        self.med_diff_data_path = None

        # Death data window and path
        self.death_data_win = None
        self.death_data_path = None

        # Controls for population data
        self.plot_pop_data = plot_pop_data_freq > 0
        self.write_pop_data = write_pop_data_freq > 0

        # Controls for difusion data
        self.plot_med_diff_data = plot_med_diff_data_freq > 0
        self.write_med_diff_data = write_med_diff_data_freq > 0

        # Controls for death data
        self.plot_death_data = plot_death_data_freq > 0
        self.write_death_data = write_death_data_freq > 0

        # Cell death mechanism tracking
        self.death_mech = {'viral': 0,
                           'contact': 0}

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
        model_string = immune_recruitment_model_string(num_ec=num_ecs,
                                                       sites_per_cell=cell_volume,
                                                       time_conv=s_to_mcs / 60 / 60 / 24)
        self.add_free_floating_antimony(model_string=model_string,
                                        model_name='TcellModel',
                                        step_size=sbml_step_size)

        #   Get reference to solver
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == 'TcellModel':
                self.rr = rr

        # Step: load secretion values from ODEs into XML

        #   Extracellular virus: diffusion and decay
        self.get_xml_element('virus_dc').cdata = virus_dc
        self.get_xml_element('virus_decay').cdata = self.rr["c"]

        #   Local cytokine: diffusion and decay
        self.get_xml_element('cytokine_dc').cdata = cytokine_dc
        self.get_xml_element('cytokine_decay').cdata = self.rr["cc"] + self.rr["kc"]

        # Step: data tracking setup

        # Initialize windows and paths if requested

        dot_size = 5
        line_size = 3

        if self.plot_pop_data:
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

            if plot_ode_sol:
                self.pop_data_win.add_plot("UninfectedODE", style='Lines', color='blue', size=line_size)
                self.pop_data_win.add_plot("InfectedODE", style='Lines', color='red', size=line_size)
                self.pop_data_win.add_plot("VirusReleasingODE", style='Lines', color='green', size=line_size)
                self.pop_data_win.add_plot("DeadODE", style='Lines', color='yellow', size=line_size)
                self.pop_data_win.add_plot("CD8LocalODE", style='Lines', color='white', size=line_size)
                self.pop_data_win.add_plot("CD8LymphODE", style='Lines', color='purple', size=line_size)

        if self.plot_med_diff_data:
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

            if plot_ode_sol:
                self.med_diff_data_win.add_plot("MedViralODE", style='Lines', color='red', size=line_size)
                self.med_diff_data_win.add_plot("MedCytLocalODE", style='Lines', color='blue', size=line_size)
                self.med_diff_data_win.add_plot("MedCytLymphODE", style='Lines', color='green', size=line_size)

        if self.plot_death_data:
            self.death_data_win = self.add_new_plot_window(title='Death data',
                                                           x_axis_title='MCS',
                                                           y_axis_title='Numer of cells',
                                                           x_scale_type='linear',
                                                           y_scale_type='log',
                                                           grid=True,
                                                           config_options={'legend': True})

            self.death_data_win.add_plot("Viral", style='Dots', color='blue', size=5)
            self.death_data_win.add_plot("Contact", style='Dots', color='green', size=5)

            if plot_ode_sol:
                self.death_data_win.add_plot("ViralODE", style='Lines', color='blue', size=line_size)
                self.death_data_win.add_plot("ContactODE", style='Lines', color='green', size=line_size)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_pop_data:
                self.pop_data_path = Path(self.output_dir).joinpath('pop_data.dat')
                with open(self.pop_data_path, 'w'):
                    pass

            if self.write_med_diff_data:
                self.med_diff_data_path = Path(self.output_dir).joinpath('med_diff_data.dat')
                with open(self.med_diff_data_path, 'w'):
                    pass

            if self.write_death_data:
                self.death_data_path = Path(self.output_dir).joinpath('death_data.dat')
                with open(self.death_data_path, 'w'):
                    pass

        if plot_ode_sol:
            #   Generate solver instance
            model_string = immune_recruitment_model_string_ode(num_ec=num_ecs,
                                                               time_conv=s_to_mcs / 60 / 60 / 24)
            self.add_free_floating_antimony(model_string=model_string,
                                            model_name='TcellModelODE',
                                            step_size=sbml_step_size)

            #   Get reference to solver
            from cc3d.CompuCellSetup import persistent_globals as pg
            for model_name, rr in pg.free_floating_sbml_simulators.items():
                if model_name == 'TcellModelODE':
                    self.rr_ode = rr

            #   Apply initial conditions
            self.rr_ode["T"] = len(self.cell_list_by_type(self.UNINFECTED))
            self.rr_ode["I1"] = len(self.cell_list_by_type(self.INFECTED))
            self.rr_ode["I2"] = len(self.cell_list_by_type(self.VIRUSRELEASING))
            self.rr_ode["D"] = len(self.cell_list_by_type(self.DEAD))

    def step(self, mcs):
        # Prep stuff
        num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))
        virus_secretor = self.get_field_secretor("Virus")
        cytokine_secretor = self.get_field_secretor("cytokine")

        # Step: Immune cell killing
        dei2 = self.rr["dei2"]

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
            pr_death = death_rate * math.exp(-death_rate)
            if random.random() < pr_death:
                cell.type = self.DEAD
                cell.targetVolume = 0
                self.death_mech['contact'] += 1

        # Step: Viral killing

        death_rate = self.rr["d"]
        pr_death = death_rate * math.exp(-death_rate)
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            if random.random() < pr_death:
                cell.type = self.DEAD
                self.death_mech['viral'] += 1

        # Step: Viral life cycle

        # Infected -> Virus releasing
        release_rate = self.rr["k"]
        pr_release = release_rate * math.exp(-release_rate)
        for cell in self.cell_list_by_type(self.INFECTED):
            if random.random() < pr_release:
                cell.type = self.VIRUSRELEASING

        # Internalization
        beta = self.rr["beta"]
        for cell in self.cell_list_by_type(self.UNINFECTED):
            seen_amount = virus_secretor.amountSeenByCell(cell) / cell.volume * self.dim.z
            infect_rate = beta * seen_amount
            pr_infect = infect_rate * math.exp(-infect_rate)
            if random.random() < pr_infect:
                cell.type = self.INFECTED

        # Step: Immune recruitment

        # Pass spatial information to ODEs
        self.rr["C"] = cytokine_secretor.totalFieldIntegral() / self.dim.z

        # Integrate ODEs
        self.rr.timestep()

        # Calculate outflow
        remove_rate = self.rr["dE"]
        pr_remove = remove_rate * math.exp(-remove_rate)
        for cell in self.cell_list_by_type(self.CD8LOCAL):
            if random.random() < pr_remove:
                cell.targetVolume = 0

        # Calculate inflow
        num_add = 0
        add_rate = self.rr["El"] * self.rr["ke"]
        exp_term = math.exp(-add_rate)
        sum_term = 1.0
        while random.random() < 1 - exp_term * sum_term:
            num_add += 1.0
            sum_term += add_rate ** num_add / math.factorial(num_add)
        print(f'self.rr["El"], self.rr["ke"], add_rate, num_add: ', self.rr["El"], self.rr["ke"], add_rate, num_add)

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
        secr_amount = self.rr["p"]
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            virus_secretor.secreteInsideCellTotalCount(cell, secr_amount / cell.volume)

        #   Local cytokine: secretion by infected and virus-releasing cells
        secr_amount = self.rr["pc"]
        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING):
            cytokine_secretor.secreteInsideCellTotalCount(cell, secr_amount / cell.volume)

        # Step: Update data tracking

        if plot_ode_sol:
            self.rr_ode.timestep()

        # Check if plotting/writing
        plot_pop_data = self.plot_pop_data and mcs % plot_pop_data_freq == 0
        plot_med_diff_data = self.plot_med_diff_data and mcs % plot_med_diff_data_freq == 0
        plot_death_data = self.plot_death_data and mcs % plot_death_data_freq == 0

        # Check for existence of output directory and disable writing if necessary
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % write_pop_data_freq == 0
            write_med_diff_data = self.write_med_diff_data and mcs % write_med_diff_data_freq == 0
            write_death_data = self.write_death_data and mcs % write_death_data_freq == 0
        else:
            write_pop_data = False
            write_med_diff_data = False
            write_death_data = False

        # Population data tracking

        if plot_pop_data or write_pop_data:

            # Gather population data
            num_cells_uninfected = len(self.cell_list_by_type(self.UNINFECTED))
            num_cells_infected = len(self.cell_list_by_type(self.INFECTED))
            num_cells_virusreleasing = len(self.cell_list_by_type(self.VIRUSRELEASING))
            num_cells_dead = len(self.cell_list_by_type(self.DEAD))
            num_cells_immune = len(self.cell_list_by_type(self.CD8LOCAL))
            num_cells_immune_l = self.rr_ode["El"]

            # Plot population data plot if requested
            min_thresh = 0.1
            if plot_pop_data:
                if num_cells_uninfected > min_thresh:
                    self.pop_data_win.add_data_point('Uninfected', mcs, num_cells_uninfected)
                if num_cells_infected > min_thresh:
                    self.pop_data_win.add_data_point('Infected', mcs, num_cells_infected)
                if num_cells_virusreleasing > min_thresh:
                    self.pop_data_win.add_data_point('VirusReleasing', mcs, num_cells_virusreleasing)
                if num_cells_dead > min_thresh:
                    self.pop_data_win.add_data_point('Dead', mcs, num_cells_dead)
                if num_cells_immune > min_thresh:
                    self.pop_data_win.add_data_point('CD8Local', mcs, num_cells_immune)
                if num_cells_immune_l > min_thresh:
                    self.pop_data_win.add_data_point('CD8Lymph', mcs, num_cells_immune_l)

                if plot_ode_sol:
                    num_cells_uninfected = self.rr_ode["T"]
                    num_cells_infected = self.rr_ode["I1"]
                    num_cells_virusreleasing = self.rr_ode["I2"]
                    num_cells_dead = self.rr_ode["D"]
                    num_cells_immune = self.rr_ode["E"]
                    num_cells_immune_l = self.rr_ode["El"]

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

            # Write population data to file if requested
            if write_pop_data:
                with open(self.pop_data_path, 'a') as fout:
                    fout.write('{}, {}, {}, {}, {}, {}\n'.format(mcs,
                                                                 num_cells_uninfected,
                                                                 num_cells_infected,
                                                                 num_cells_virusreleasing,
                                                                 num_cells_dead,
                                                                 num_cells_immune))

        # Diffusive field data tracking

        if plot_med_diff_data or write_med_diff_data:

            # Gather total diffusive amounts
            med_viral_total = virus_secretor.totalFieldIntegral() / self.dim.z
            med_cyt_total = cytokine_secretor.totalFieldIntegral() / self.dim.z
            med_cyt_lymph = self.rr["Cl"]
            # Plot total diffusive viral amount if requested
            if plot_med_diff_data:
                if med_viral_total > 0:
                    self.med_diff_data_win.add_data_point("MedViral", mcs, med_viral_total)
                if med_cyt_total > 0:
                    self.med_diff_data_win.add_data_point("MedCytLocal", mcs, med_cyt_total)
                if med_cyt_lymph > 0:
                    self.med_diff_data_win.add_data_point("MedCytLymph", mcs, med_cyt_lymph)

            # Write total diffusive viral amount if requested
            if write_med_diff_data:
                with open(self.med_diff_data_path, 'a') as fout:
                    fout.write('{}, {}, {}, {}\n'.format(mcs, med_viral_total, med_cyt_total, med_cyt_lymph))

            if plot_ode_sol:
                med_viral_total = self.rr_ode["V"]
                med_cyt_total = self.rr_ode["C"]
                med_cyt_lymph = self.rr_ode["Cl"]

                # Plot total diffusive viral amount if requested
                if plot_med_diff_data:
                    if med_viral_total > 0:
                        self.med_diff_data_win.add_data_point("MedViralODE", mcs, med_viral_total)
                    if med_cyt_total > 0:
                        self.med_diff_data_win.add_data_point("MedCytLocalODE", mcs, med_cyt_total)
                    if med_cyt_lymph > 0:
                        self.med_diff_data_win.add_data_point("MedCytLymphODE", mcs, med_cyt_lymph)

        # Death mechanism data tracking

        if plot_death_data or write_death_data:
            num_viral = self.death_mech['viral']
            num_contact = self.death_mech['contact']

            # Plot death data if requested
            min_thresh = 0.1
            if plot_death_data:
                if num_viral > min_thresh:
                    self.death_data_win.add_data_point("Viral", mcs, num_viral)
                if num_contact > min_thresh:
                    self.death_data_win.add_data_point("Contact", mcs, num_contact)

            # Write death data if requested
            if write_death_data:
                with open(self.death_data_path, 'a') as fout:
                    fout.write(f'{mcs}, {num_viral}, {num_contact}\n')

            if plot_ode_sol:
                num_viral = self.rr_ode['viralDeath']
                num_contact = self.rr_ode['cd8Death']

                # Plot death data if requested
                if plot_death_data:
                    if num_viral > min_thresh:
                        self.death_data_win.add_data_point("ViralODE", mcs, num_viral)
                    if num_contact > min_thresh:
                        self.death_data_win.add_data_point("ContactODE", mcs, num_contact)

