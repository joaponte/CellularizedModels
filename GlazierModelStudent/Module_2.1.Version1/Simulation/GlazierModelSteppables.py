from cc3d.core.PySteppables import *

import math
import random

# Conversion Factors
s_to_mcs = 5.0 * 60.0  # s/mcs
days_to_mcs = s_to_mcs / 60.0 / 60.0 / 24.0  # days/mcs

# =============================
# CompuCell Parameters
# cell
cell_diameter = 5.0
cell_volume = cell_diameter ** 2
volume_lm = 9

# Antimony/SBML model step size
sbml_step_size = days_to_mcs 

ode_model_string = '''model ODEModel()
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

class GlazierModelSteppable(SteppableBasePy):

    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        
        # Population data window and path
        self.pop_data_win = None 
            
        # Diffusion data window and path
        self.med_diff_data_win = None
        
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
        
        # Add floating antimony model of the ODE model
        self.add_free_floating_antimony(model_string=ode_model_string,
                                        model_name='ODEModel',
                                        step_size=sbml_step_size)
                                        
        # Step: calculate scaling coefficients
        num_ecs = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DEAD))
        self.scale_by_volume = num_ecs / self.sbml.ODEModel["T0"]
        self.scale_by_voxel = 1 / self.sbml.ODEModel["T0"] / cell_volume
        self.scale_by_time = days_to_mcs
        
        # Population data tracking plot setup
        self.pop_data_win = self.add_new_plot_window(title='Population data',
                                                     x_axis_title='Time (days)',
                                                     y_axis_title='Numer of cells',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
        
        dot_size = 5
        self.pop_data_win.add_plot("Uninfected", style='Dots', color='blue', size=dot_size)
        self.pop_data_win.add_plot("Infected", style='Dots', color='red', size=dot_size)
        self.pop_data_win.add_plot("VirusReleasing", style='Dots', color='green', size=dot_size)
        self.pop_data_win.add_plot("Dead", style='Dots', color='yellow', size=dot_size)
        
        line_size = 3
        self.pop_data_win.add_plot("UninfectedODE", style='Lines', color='blue', size=line_size)
        self.pop_data_win.add_plot("InfectedODE", style='Lines', color='red', size=line_size)
        self.pop_data_win.add_plot("VirusReleasingODE", style='Lines', color='green', size=line_size)
        self.pop_data_win.add_plot("DeadODE", style='Lines', color='yellow', size=line_size)
        self.pop_data_win.add_plot("CD8LocalODE", style='Lines', color='white', size=line_size)
        self.pop_data_win.add_plot("CD8LymphODE", style='Lines', color='purple', size=line_size)
        
        # Diffusive field data tracking setup
        self.med_diff_data_win = self.add_new_plot_window(title='Total diffusive species',
                                                          x_axis_title='Time (days)',
                                                          y_axis_title='Number of diffusive species per volume',
                                                          x_scale_type='linear',
                                                          y_scale_type='log',
                                                          grid=True,
                                                          config_options={'legend': True})
        
        self.med_diff_data_win.add_plot("MedViralODE", style='Lines', color='red', size=line_size)
        self.med_diff_data_win.add_plot("MedCytLocalODE", style='Lines', color='blue', size=line_size)
        self.med_diff_data_win.add_plot("MedCytLymphODE", style='Lines', color='green', size=line_size)
        
        # Death mechanism data tracking setup
        self.death_data_win = self.add_new_plot_window(title='Death data',
                                                       x_axis_title='Time (days)',
                                                       y_axis_title='Numer of cells',
                                                       x_scale_type='linear',
                                                       y_scale_type='linear',
                                                       grid=True,
                                                       config_options={'legend': True})
                                                       
        
        self.death_data_win.add_plot("Viral", style='Dots', color='blue', size=5)
        
        self.death_data_win.add_plot("ViralODE", style='Lines', color='blue', size=line_size)
        self.death_data_win.add_plot("ContactODE", style='Lines', color='green', size=line_size)
        
    def step(self,mcs):
        # Step: Viral killing
        death_rate = self.sbml.ODEModel["d"] * self.scale_by_time
        pr_death = 1 - math.exp(-death_rate)
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            if random.random() < pr_death:
                cell.type = self.DEAD
                self.death_mech['viral'] += 1
        
        num_cells_dead = len(self.cell_list_by_type(self.DEAD))
        self.pop_data_win.add_data_point('Dead', mcs * days_to_mcs, num_cells_dead)
        
        # Plot viral death data
        if self.death_mech['viral'] > 0:
            self.death_data_win.add_data_point("Viral", mcs * days_to_mcs, self.death_mech['viral'])
        
        # Infected -> Virus releasing
        release_rate = self.sbml.ODEModel["k"] * self.scale_by_time
        pr_release = 1 - math.exp(-release_rate)
        for cell in self.cell_list_by_type(self.INFECTED):
            if random.random() < pr_release:
                cell.type = self.VIRUSRELEASING
                
        num_cells_virusreleasing = len(self.cell_list_by_type(self.VIRUSRELEASING))
        self.pop_data_win.add_data_point('VirusReleasing', mcs * days_to_mcs, num_cells_virusreleasing)
        
        # Uninfected -> Infected
        beta = self.sbml.ODEModel["beta"] * self.scale_by_time
        infect_rate = beta * self.sbml.ODEModel["V"]
        pr_infect = 1 - math.exp(-infect_rate)
        for cell in self.cell_list_by_type(self.UNINFECTED):
            if random.random() < pr_infect:
                cell.type = self.INFECTED
        
        num_cells_uninfected = len(self.cell_list_by_type(self.UNINFECTED))
        self.pop_data_win.add_data_point('Uninfected', mcs * days_to_mcs, num_cells_uninfected)
        num_cells_infected = len(self.cell_list_by_type(self.INFECTED))
        self.pop_data_win.add_data_point('Infected', mcs * days_to_mcs, num_cells_infected)
        
        # Integrate ODEs
        self.timestep_sbml()

        # Population data tracking
        num_cells_uninfected = self.sbml.ODEModel["T"] * self.scale_by_volume
        num_cells_infected = self.sbml.ODEModel["I1"] * self.scale_by_volume
        num_cells_virusreleasing = self.sbml.ODEModel["I2"] * self.scale_by_volume
        num_cells_dead = self.sbml.ODEModel["D"] * self.scale_by_volume
        num_cells_immune = self.sbml.ODEModel["E"] * self.scale_by_volume
        num_cells_immune_l = self.sbml.ODEModel["El"] * self.scale_by_volume

        self.pop_data_win.add_data_point('UninfectedODE', mcs * days_to_mcs, num_cells_uninfected)
        self.pop_data_win.add_data_point('InfectedODE', mcs * days_to_mcs, num_cells_infected)
        self.pop_data_win.add_data_point('VirusReleasingODE', mcs * days_to_mcs, num_cells_virusreleasing)
        self.pop_data_win.add_data_point('DeadODE', mcs * days_to_mcs, num_cells_dead)
#         self.pop_data_win.add_data_point('CD8LocalODE', mcs * days_to_mcs, num_cells_immune)
#         self.pop_data_win.add_data_point('CD8LymphODE', mcs * days_to_mcs, num_cells_immune_l)
        
        # Diffusive field data tracking
        med_viral_total = self.sbml.ODEModel["V"] * self.scale_by_volume
        med_cyt_total = self.sbml.ODEModel["C"] * self.scale_by_volume
        med_cyt_lymph = self.sbml.ODEModel["Cl"] * self.scale_by_volume
        
        if med_viral_total > 0:
            self.med_diff_data_win.add_data_point("MedViralODE", mcs * days_to_mcs, med_viral_total)
        if med_cyt_total > 0:
            self.med_diff_data_win.add_data_point("MedCytLocalODE", mcs * days_to_mcs, med_cyt_total)
        if med_cyt_lymph > 0:
            self.med_diff_data_win.add_data_point("MedCytLymphODE", mcs * days_to_mcs, med_cyt_lymph)

        # Death mechanism data tracking
        num_viral = self.sbml.ODEModel['viralDeath'] * self.scale_by_volume
        num_contact = self.sbml.ODEModel['cd8Death'] * self.scale_by_volume

        self.death_data_win.add_data_point("ViralODE", mcs * days_to_mcs, num_viral)
        self.death_data_win.add_data_point("ContactODE", mcs * days_to_mcs, num_contact)

    def finish(self):
        """
        Finish Function is called after the last MCS
        """


        