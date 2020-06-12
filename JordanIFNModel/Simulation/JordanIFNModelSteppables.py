from cc3d.core.PySteppables import *
import numpy as np
import os

min_to_mcs = 10.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0  # hrs/mcs
hours_to_simulate = 36.0

plot_ODEModel = True
plot_CellularizedModel = True

## How_to_determine_IFNe
# 1 determines IFNe from ODE model
# 2 determines IFNe from scalar from CC3D
# 3 determines IFNe from field from CC3D
How_to_determine_IFNe = 1

'''Jordan J. A. Weaver and Jason E. Shoemaker. Mathematical Modeling of RNA Virus Sensing Pathways Reveal Paracrine Signaling as the Primary Factor 
Regulating Excessive Cytokine Production'''

# Modified code to make IFNe production dependent on the number of infected cells (P)
model_string = '''
    //Equations
    E2a: -> IFN         ; P*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)    ;
    E2b: IFN -> IFNe    ; k21*IFN                                           ;
    E3a: IFNe ->        ; t2*IFNe                                           ;
    E4a: -> STATP       ; P*k31*IFNe/(k32+k33*IFNe)                         ;
    E4b: STATP ->       ; t3*STATP                                          ;
    E5a: -> IRF7        ; P*(k41*STATP+k42*IRF7P)                           ;
    E5b: IRF7 ->        ; t4*IRF7                                           ;
    E6a: -> IRF7P       ; P*k51*IRF7                                        ;
    E6b: IRF7P ->       ; t5*IRF7P                                          ;
    E7a: P ->           ; P*k61*V                                           ;
    E8a: -> V           ; P*(k71*V)/(1.0+k72*IFN*7E-5)                      ;
    E8b: V ->           ; k73*V                                             ;
    
    //Parameters
    k11 = 0.0       ; 
    k12 = 9.746     ; 
    k13 = 12.511    ; 
    k14 = 13.562    ;
    k21 = 10.385    ;
    t2  = 3.481     ;
    k31 = 45.922    ;
    k32 = 5.464     ;
    k33 = 0.068     ;
    t3  = 0.3       ;
    k41 = 0.115     ;
    k42 = 1.053     ;
    t4  = 0.3       ;
    k51 = 0.202     ;
    t5  = 0.3       ;
    k61 = 0.635     ;
    k71 = 1.537     ;
    k72 = 47.883    ;
    k73 = 0.197     ;
    n   = 3.0       ;

    //Initial Conditions
    P    = 1.0      ;
    RIGI = 1.0      ;
    IRF7 = 0.72205  ;
    V    = 6.9e-8   ;
'''
class ODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Uptading max simulation steps using scaling factor to simulate 10 days
        self.get_xml_element('simulation_steps').cdata = hours_to_simulate / hours_to_mcs

        # Adding free floating antimony model
        self.add_free_floating_antimony(model_string=model_string, model_name='ODEModel',
                                        step_size=hours_to_mcs)
    
    def step(self,mcs):
        self.timestep_sbml()

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # set initial model parameters
        self.get_xml_element('IFNe_decay').cdata = self.sbml.ODEModel['t2'] * hours_to_mcs
        self.initial_infected = len(self.cell_list_by_type(self.I2))
        self.shared_steppable_vars['IFNe'] = self.sbml.ODEModel['IFNe']
        self.shared_steppable_vars['STATP'] = self.sbml.ODEModel['STATP']
        self.shared_steppable_vars['IRF7'] = self.sbml.ODEModel['IRF7']
        self.shared_steppable_vars['IRF7P'] = self.sbml.ODEModel['IRF7P']
        self.shared_steppable_vars['V'] = self.sbml.ODEModel['V']

    def step(self, mcs):
        pIFNe = 0.0
        pSTATP = 0.0
        pIRF7 = 0.0
        pIRF7P = 0.0
        pV = 0.0
        secretor = self.get_field_secretor("IFNe")
        for cell in self.cell_list_by_type(self.I2):
            # Rule 3a
            k21 = self.sbml.ODEModel['k21'] * hours_to_mcs / self.initial_infected
            IFN = self.sbml.ODEModel['IFN']
            release = secretor.secreteInsideCellTotalCount(cell, k21 * IFN / cell.volume)
            pIFNe += abs(release.tot_amount)

            # Rule 4a
            k31 = self.sbml.ODEModel['k31'] * hours_to_mcs / self.initial_infected
            IFNe = self.shared_steppable_vars['IFNe']
            k32 = self.sbml.ODEModel['k32']
            k33 = self.sbml.ODEModel['k33']
            pSTATP += k31*IFNe/(k32+k33*IFNe)

            # Rule 5a
            k41 = self.sbml.ODEModel['k41'] * hours_to_mcs / self.initial_infected
            STATP = self.shared_steppable_vars['STATP']
            k42 = self.sbml.ODEModel['k42'] * hours_to_mcs / self.initial_infected
            IRF7P = self.shared_steppable_vars['IRF7P']
            pIRF7 += k41 * STATP + k42 * IRF7P

            # Rule 6a
            k51 = self.sbml.ODEModel['k51'] * hours_to_mcs / self.initial_infected
            IRF7 = self.shared_steppable_vars['IRF7']
            pIRF7P += k51 * IRF7

            # Rule 7a
            k61 = self.sbml.ODEModel['k61'] * hours_to_mcs
            V = self.sbml.ODEModel['V']
            p_I2toDead = k61 * V
            if np.random.random() < p_I2toDead:
                cell.type = self.DEAD

            # Rule 8a
            k71 = self.sbml.ODEModel['k71'] * hours_to_mcs / self.initial_infected * 1.1 # DONT UNDERSTAND WHY
            V = self.shared_steppable_vars['V']
            k72 = self.sbml.ODEModel['k72']
            pV += k71 * V / (1.0 + k72 * IFN * 7E-5)

        # Rule 3b
        t2 = self.sbml.ODEModel['t2'] * hours_to_mcs
        self.shared_steppable_vars['IFNe'] += pIFNe - t2 * self.shared_steppable_vars['IFNe']

        # Rule 4b
        t3 = self.sbml.ODEModel['t3'] * hours_to_mcs
        self.shared_steppable_vars['STATP'] += pSTATP - t3 * self.shared_steppable_vars['STATP']

        # Rule 5b
        t4 = self.sbml.ODEModel['t4'] * hours_to_mcs
        self.shared_steppable_vars['IRF7'] += pIRF7 - t4 * self.shared_steppable_vars['IRF7']

        # Rule 6b
        t5 = self.sbml.ODEModel['t5'] * hours_to_mcs
        self.shared_steppable_vars['IRF7P'] += pIRF7P - t5 * self.shared_steppable_vars['IRF7P']

        # Rule 8b
        k73 = self.sbml.ODEModel['k73'] * hours_to_mcs
        self.shared_steppable_vars['V'] += pV - k73 * self.shared_steppable_vars['V']

class PlotODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_infected = len(self.cell_list_by_type(self.I2))

        # Initialize Graphic Window for Amber Smith ODE model
        if (plot_ODEModel or plot_CellularizedModel):
            self.plot_win = self.add_new_plot_window(title='P',
                                                     x_axis_title='Hours',
                                                     y_axis_title='Number of Cells', x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False, config_options={'legend': True})

            self.plot_win2 = self.add_new_plot_window(title='STATP',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win3 = self.add_new_plot_window(title='IFNe',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Extracellular IFN', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})


            self.plot_win4 = self.add_new_plot_window(title='IRF7',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win5 = self.add_new_plot_window(title='IRF7P',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win6 = self.add_new_plot_window(title='V',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            if plot_ODEModel:
                self.plot_win.add_plot("ODEP", style='Dots', color='red', size=5)
                self.plot_win2.add_plot("ODESTATP", style='Dots', color='blue', size=5)
                self.plot_win3.add_plot("ODEIFNe", style='Dots', color='purple', size=5)
                self.plot_win4.add_plot("ODEIRF7", style='Dots', color='orange', size=5)
                self.plot_win5.add_plot("ODEIRF7P", style='Dots', color='green', size=5)
                self.plot_win6.add_plot("ODEV", style='Dots', color='yellow', size=5)

            if plot_CellularizedModel:
                self.plot_win.add_plot("CC3DP", style='Lines', color='red', size=5)
                self.plot_win2.add_plot("CC3DSTATP", style='Lines', color='blue', size=5)
                self.plot_win3.add_plot("CC3DIFNe", style='Lines', color='purple', size=5)
                self.plot_win4.add_plot("CC3DIRF7", style='Lines', color='orange', size=5)
                self.plot_win5.add_plot("CC3DIRF7P", style='Lines', color='green', size=5)
                self.plot_win6.add_plot("CC3DV", style='Lines', color='yellow', size=5)

    def step(self, mcs):
        if plot_ODEModel:
            self.plot_win.add_data_point("ODEP", mcs * hours_to_mcs,self.sbml.ODEModel['P'])
            self.plot_win2.add_data_point("ODESTATP", mcs * hours_to_mcs, self.sbml.ODEModel['STATP'])
            self.plot_win3.add_data_point("ODEIFNe", mcs * hours_to_mcs, self.sbml.ODEModel['IFNe'])
            self.plot_win4.add_data_point("ODEIRF7", mcs * hours_to_mcs, self.sbml.ODEModel['IRF7'])
            self.plot_win5.add_data_point("ODEIRF7P", mcs * hours_to_mcs, self.sbml.ODEModel['IRF7P'])
            self.plot_win6.add_data_point("ODEV", mcs * hours_to_mcs, self.sbml.ODEModel['V'])

        if plot_CellularizedModel:
            P = len(self.cell_list_by_type(self.I2))/self.initial_infected
            self.plot_win.add_data_point("CC3DP", mcs * hours_to_mcs, P)
            self.plot_win2.add_data_point("CC3DSTATP", mcs * hours_to_mcs,
                                          self.shared_steppable_vars['STATP'])
            self.plot_win3.add_data_point("CC3DIFNe", mcs * hours_to_mcs,
                                          self.shared_steppable_vars['IFNe'])
            self.plot_win4.add_data_point("CC3DIRF7", mcs * hours_to_mcs,
                                          self.shared_steppable_vars['IRF7'])
            self.plot_win5.add_data_point("CC3DIRF7P", mcs * hours_to_mcs,
                                          self.shared_steppable_vars['IRF7P'])
            self.plot_win6.add_data_point("CC3DV", mcs * hours_to_mcs,
                                          self.shared_steppable_vars['V'])