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
        
    def start(self):
        self.add_free_floating_antimony(model_string=ode_model_string,
                                        model_name='ODEModel',
                                        step_size=sbml_step_size)
                                        
        
        # Population data tracking plot setup
        self.pop_data_win = self.add_new_plot_window(title='Population data',
                                                     x_axis_title='Time (days)',
                                                     y_axis_title='Numer of cells',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})
        line_size = 3
        self.pop_data_win.add_plot("UninfectedODE", style='Lines', color='blue', size=line_size)
        self.pop_data_win.add_plot("InfectedODE", style='Lines', color='red', size=line_size)
        self.pop_data_win.add_plot("VirusReleasingODE", style='Lines', color='green', size=line_size)
        self.pop_data_win.add_plot("DeadODE", style='Lines', color='yellow', size=line_size)
        self.pop_data_win.add_plot("CD8LocalODE", style='Lines', color='white', size=line_size)
        self.pop_data_win.add_plot("CD8LymphODE", style='Lines', color='purple', size=line_size)
        
        
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

    def step(self,mcs):
        # Integrate ODEs
        self.timestep_sbml()
        
        # Population data tracking
        num_cells_uninfected = self.sbml.ODEModel["T"]
        num_cells_infected = self.sbml.ODEModel["I1"]
        num_cells_virusreleasing = self.sbml.ODEModel["I2"]
        num_cells_dead = self.sbml.ODEModel["D"]
        num_cells_immune = self.sbml.ODEModel["E"]
        num_cells_immune_l = self.sbml.ODEModel["El"]

        self.pop_data_win.add_data_point('UninfectedODE', mcs * days_to_mcs, num_cells_uninfected)
        self.pop_data_win.add_data_point('InfectedODE', mcs * days_to_mcs, num_cells_infected)
        self.pop_data_win.add_data_point('VirusReleasingODE', mcs * days_to_mcs, num_cells_virusreleasing)
        self.pop_data_win.add_data_point('DeadODE', mcs * days_to_mcs, num_cells_dead)
        self.pop_data_win.add_data_point('CD8LocalODE', mcs * days_to_mcs, num_cells_immune)
        self.pop_data_win.add_data_point('CD8LymphODE', mcs * days_to_mcs, num_cells_immune_l)
        
        # Diffusive field data tracking
        med_viral_total = self.sbml.ODEModel["V"]
        med_cyt_total = self.sbml.ODEModel["C"]
        med_cyt_lymph = self.sbml.ODEModel["Cl"]
        
        if med_viral_total > 0:
            self.med_diff_data_win.add_data_point("MedViralODE", mcs * days_to_mcs, med_viral_total)
        if med_cyt_total > 0:
            self.med_diff_data_win.add_data_point("MedCytLocalODE", mcs * days_to_mcs, med_cyt_total)
        if med_cyt_lymph > 0:
            self.med_diff_data_win.add_data_point("MedCytLymphODE", mcs * days_to_mcs, med_cyt_lymph)

    def finish(self):
        """
        Finish Function is called after the last MCS
        """


        