from cc3d.core.PySteppables import *

import os
import sys

import math
import random

# Import project libraries and classes
sys.path.append(os.path.dirname(__file__))
import GlazierModelLib
from GlazierModelInputs import *


class CellsInitializerSteppable(SteppableBasePy):
    """
    Initializes epithelial sheet and an initial immune cell population
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
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
        if init_infect == 0:
            cell = self.cell_field[self.dim.x // 2, self.dim.y // 2, 0]
            cell.type = self.INFECTED
        else:
            num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))
            num_to_infect = int(init_infect * num_ecs)
            num_infected = 0
            while num_infected < num_to_infect:
                xi = random.randint(0, self.dim.x)
                yi = random.randint(0, self.dim.y)
                cell = self.cell_field[xi, yi, 0]
                if cell is not None and cell.type == self.UNINFECTED:
                    cell.type = self.INFECTED
                    num_infected += 1


class ViralSecretionSteppable(SteppableBasePy):
    """
    Implements viral release module
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to ImmuneRecruitmentSteppable
        self.ir_steppable = None

    def start(self):
        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[GlazierModelLib.ir_steppable_key]

        self.get_xml_element('virus_dc').cdata = virus_dc
        self.get_xml_element('virus_decay').cdata = self.ir_steppable.get_model_val("c")

        secr_amount = self.ir_steppable.get_model_val('p') / cell_volume
        self.get_xml_element("virus_secr3").cdata = secr_amount


class ImmuneCellKillingSteppable(SteppableBasePy):
    """
    Implements immune cell direct cytotoxicity and bystander effect module
    """
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        # Reference to ImmuneRecruitmentSteppable
        self.ir_steppable = None

        self.num_ecs = None

    def start(self):
        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[GlazierModelLib.ir_steppable_key]

        self.num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))

    def step(self, mcs):
        if self.simdata_steppable is None:
            self.simdata_steppable: SimDataSteppable = \
                self.shared_steppable_vars[GlazierModelLib.simdata_steppable_key]

        dei2 = self.ir_steppable.get_model_val("dei2")

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

            death_rate = dei2 * self.num_ecs / srf_area * cd8_area
            pr_death = GlazierModelLib.ul_rate_to_prob(death_rate)
            if random.random() < pr_death:
                cell.type = self.DEAD
                cell.targetVolume = 0
                self.simdata_steppable.track_death_contact()


class ChemotaxisSteppable(SteppableBasePy):
    """
    Implements immune cell chemotaxis module
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        for cell in self.cell_list_by_type(self.CD8LOCAL):
            cd = self.chemotaxisPlugin.addChemotaxisData(cell, "cytokine")
            cd.setLambda(lamda_chemotaxis)
            cd.assignChemotactTowardsVectorTypes([self.MEDIUM])

    def step(self, mcs):
        field = self.field.cytokine
        for cell in self.cell_list_by_type(self.CD8LOCAL):
            cd = self.chemotaxisPlugin.getChemotaxisData(cell, "cytokine")
            if cd is None:
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "cytokine")
                cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
            concentration = field[cell.xCOM, cell.yCOM, 1]
            cd.setLambda(lamda_chemotaxis / (1.0 + concentration))


class ViralLifeCycleSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        # Reference to ImmuneRecruitmentSteppable
        self.ir_steppable = None

    def start(self):

        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[GlazierModelLib.ir_steppable_key]

    def step(self, mcs):
        if self.simdata_steppable is None:
            self.simdata_steppable: SimDataSteppable = \
                self.shared_steppable_vars[GlazierModelLib.simdata_steppable_key]

        secretor = self.get_field_secretor("Virus")

        # Virus releasing -> Dead
        d = self.ir_steppable.get_model_val("d")
        pr_death = GlazierModelLib.ul_rate_to_prob(d)
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            if random.random() < pr_death:
                cell.type = self.DEAD
                self.simdata_steppable.track_death_viral()

        # Infected -> Virus releasing
        k = self.ir_steppable.get_model_val("k")
        pr_release = GlazierModelLib.ul_rate_to_prob(k)
        for cell in self.cell_list_by_type(self.INFECTED):
            if random.random() < pr_release:
                cell.type = self.VIRUSRELEASING

        # Internalization
        beta = self.ir_steppable.get_model_val("beta")
        for cell in self.cell_list_by_type(self.UNINFECTED):
            seen_amount = secretor.amountSeenByCell(cell) / cell.volume * self.dim.z
            pr_infect = GlazierModelLib.ul_rate_to_prob(beta * seen_amount)
            if random.random() < pr_infect:
                cell.type = self.INFECTED


class SimDataSteppable(SteppableBasePy):
    """
    Plots/writes simulation data of interest
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to ImmuneRecruitmentSteppable
        self.ir_steppable = None

        self.pop_data_win = None
        self.pop_data_path = None

        self.med_diff_data_win = None
        self.med_diff_data_path = None

        self.death_data_win = None
        self.death_data_path = None

        self.plot_pop_data = plot_pop_data_freq > 0
        self.write_pop_data = write_pop_data_freq > 0

        self.plot_med_diff_data = plot_med_diff_data_freq > 0
        self.write_med_diff_data = write_med_diff_data_freq > 0
        self.med_diff_key = "MedDiff"

        self.plot_death_data = plot_death_data_freq > 0
        self.write_death_data = write_death_data_freq > 0

        # Cell death mechanism tracking
        self.__death_mech = {'viral': 0,
                             'contact': 0}

        # ODE instance
        self.__rr = None

    def start(self):
        # Post reference to self
        self.shared_steppable_vars[GlazierModelLib.simdata_steppable_key] = self

        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[GlazierModelLib.ir_steppable_key]

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
            num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))
            model_string = GlazierModelLib.immune_recruitment_model_string_ode(num_ec=num_ecs,
                                                                               time_conv=s_to_mcs / 60 / 60 / 24)
            self.add_free_floating_antimony(model_string=model_string,
                                            model_name=GlazierModelLib.ir_model_name_ode,
                                            step_size=sbml_step_size)

            #   Get reference to solver
            from cc3d.CompuCellSetup import persistent_globals as pg
            for model_name, rr in pg.free_floating_sbml_simulators.items():
                if model_name == GlazierModelLib.ir_model_name_ode:
                    self.__rr = rr

            self.__rr["T"] = len(self.cell_list_by_type(self.UNINFECTED))
            self.__rr["I1"] = len(self.cell_list_by_type(self.INFECTED))
            self.__rr["I2"] = len(self.cell_list_by_type(self.VIRUSRELEASING))
            self.__rr["D"] = len(self.cell_list_by_type(self.DEAD))

    def step(self, mcs):

        if plot_ode_sol:
            self.__rr.timestep()

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
            num_cells_immune_l = self.ir_steppable.get_model_val("El")

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

            # Write population data to file if requested
            if write_pop_data:
                with open(self.pop_data_path, 'a') as fout:
                    fout.write('{}, {}, {}, {}, {}, {}, {}\n'.format(mcs,
                                                                     num_cells_uninfected,
                                                                     num_cells_infected,
                                                                     num_cells_virusreleasing,
                                                                     num_cells_dead,
                                                                     num_cells_immune,
                                                                     num_cells_immune_l))

            if plot_pop_data and plot_ode_sol:
                num_cells_uninfected = self.__rr["T"]
                num_cells_infected = self.__rr["I1"]
                num_cells_virusreleasing = self.__rr["I2"]
                num_cells_dead = self.__rr["D"]
                num_cells_immune = self.__rr["E"]
                num_cells_immune_l = self.__rr["El"]

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

        if plot_med_diff_data or write_med_diff_data:

            # Gather total diffusive amounts
            med_viral_total = self.get_field_secretor("Virus").totalFieldIntegral() / self.dim.z
            med_cyt_total = self.get_field_secretor("cytokine").totalFieldIntegral() / self.dim.z
            med_cyt_lymph = self.ir_steppable.get_model_val("Cl")
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

            if plot_med_diff_data and plot_ode_sol:
                med_viral_total = self.__rr["V"]
                med_cyt_total = self.__rr["C"]
                med_cyt_lymph = self.__rr["Cl"]

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
            num_viral = self.__death_mech['viral']
            num_contact = self.__death_mech['contact']

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

            if plot_death_data and plot_ode_sol:
                num_viral = self.__rr['viralDeath']
                num_contact = self.__rr['cd8Death']

                # Plot death data if requested
                if plot_death_data:
                    if num_viral > min_thresh:
                        self.death_data_win.add_data_point("ViralODE", mcs, num_viral)
                    if num_contact > min_thresh:
                        self.death_data_win.add_data_point("ContactODE", mcs, num_contact)

    def track_death_viral(self):
        self.__death_mech['viral'] += 1

    def track_death_contact(self):
        self.__death_mech['contact'] += 1


class CytokineSecretionSteppable(SteppableBasePy):
    """
    Implements cytokine production/secretion and immune cell activation module
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to ImmuneRecruitmentSteppable
        self.ir_steppable = None

    def start(self):
        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[GlazierModelLib.ir_steppable_key]

        # cytokine diff parameters
        self.get_xml_element('cytokine_dc').cdata = cytokine_dc
        cc = self.ir_steppable.get_model_val("cc")
        kc = self.ir_steppable.get_model_val("kc")
        self.get_xml_element('cytokine_decay').cdata = cc + kc

        secr_amount = self.ir_steppable.get_model_val("pc") / cell_volume
        self.get_xml_element('cytokine_secr2').cdata = secr_amount
        self.get_xml_element('cytokine_secr3').cdata = secr_amount


class ImmuneRecruitmentSteppable(SteppableBasePy):
    """
    Implements immune cell recruitment module
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        # Reference to solver
        self.__rr = None

    def start(self):
        # Post reference to self
        self.shared_steppable_vars[GlazierModelLib.ir_steppable_key] = self

        # Initialize model

        #   Generate solver instance
        num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))
        model_string = GlazierModelLib.immune_recruitment_model_string(num_ec=num_ecs,
                                                                       sites_per_cell=cell_volume,
                                                                       time_conv=s_to_mcs / 60 / 60 / 24)
        self.add_free_floating_antimony(model_string=model_string,
                                        model_name=GlazierModelLib.ir_model_name,
                                        step_size=sbml_step_size)

        #   Get reference to solver
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == GlazierModelLib.ir_model_name:
                self.__rr = rr

    def step(self, mcs):
        if mcs == 0:
            # Add initial immune cell population, if any
            E_init = self.get_model_val("E")
            if E_init is not None:
                [self.add_immune_cell() for _ in range(int(E_init))]

        # Pass spatial information to ODEs
        self.__rr["C"] = self.get_field_secretor("cytokine").totalFieldIntegral() / self.dim.z

        # Integrate ODEs
        self.__rr.timestep()

        # Calculate outflow
        remove_rate = self.get_model_val("dE")
        pr_remove = GlazierModelLib.ul_rate_to_prob(remove_rate)
        for cell in self.cell_list_by_type(self.CD8LOCAL):
            if random.random() < pr_remove:
                cell.targetVolume = 0

        # Calculate inflow
        num_add = 0
        add_rate = self.get_model_val("El") * self.get_model_val("ke")
        exp_term = math.exp(-add_rate)
        sum_term = 1.0
        while random.random() < 1 - exp_term * sum_term:
            num_add += 1.0
            sum_term += add_rate ** num_add / math.factorial(num_add)
        [self.add_immune_cell() for _ in range(int(num_add))]

    def add_immune_cell(self):
        # Accumulate unoccupied samples equal to a fraction of the number of sites and place the cell at the site
        # with the maximum value

        sample_frac = 0.01  # Move to inputs
        n_sites = self.dim.x * self.dim.y
        n_sites_frac = int(n_sites * sample_frac)

        max_concentration = 0
        target_field = self.field.cytokine
        x_seed = None
        y_seed = None
        sites_sampled = 0

        med_pixel_set = [ptd.pixel for ptd in self.pixel_tracker_plugin.getMediumPixelSet()]
        random.shuffle(med_pixel_set)

        for pixel in med_pixel_set:
            xi = pixel.x
            yi = pixel.y

            # No placing cells over a boundary
            if xi + int(cell_diameter / 2) >= self.dim.x or yi + int(cell_diameter / 2) >= self.dim.y:
                continue

            open_space = True
            for x in range(xi, xi + int(cell_diameter / 2)):
                for y in range(yi, yi + int(cell_diameter / 2)):
                    if self.cell_field[x, y, 1]:
                        open_space = False
                        break
            if open_space:
                concentration_iteration = target_field[xi, yi, 1]
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

    def get_model_val(self, _var_str):
        try:
            return self.__rr[_var_str]
        except RuntimeError:
            return None
