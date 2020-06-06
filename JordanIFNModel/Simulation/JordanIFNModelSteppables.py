from cc3d.core.PySteppables import *
import numpy as np
import os

min_to_mcs = 20.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0  # hrs/mcs
hours_to_simulate = 36.0

plot_ODEModel = True
plot_CellularizedModel = True
feedback = False

'''Jordan J. A. Weaver and Jason E. Shoemaker. Mathematical Modeling of RNA Virus Sensing Pathways Reveal Paracrine Signaling as the Primary Factor 
Regulating Excessive Cytokine Production'''

# Modified code to make IFNe production dependent on the number of infected cells (P)
model_string = '''
    //Equations
    E2a: -> IFN         ; P*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)            ; // Intracellular IFN
    E2b: IFN -> IFNe    ; k21*IFN                                                 ;
    E3: IFNe ->         ; t2*IFNe                                                   ; // Extracellular IFN
    E4: -> STATP        ; P*k31*IFNe/(k32+k33*IFNe)                                 ; // Intracellular STATP
    E4b: STATP ->       ; t3*STATP                                                  ;
    E5a: -> IRF7        ; P*(k41*STATP+k42*IRF7P)                                   ; // Intracellular IRF7
    E5b: IRF7 ->        ; t4*IRF7
    E6a: -> IRF7P       ; P*k51*IRF7                                                ; // Intracellular IRF7P
    E6b: IRF7P ->       ; t5*IRF7P
    E7a: P ->           ; P*k61*V                                                   ; // Infected Cells
    E8a: -> V           ; P*(k71*V)/(1.0+k72*IFN*7E-5)                              ; // Intracellular Virus
    E8b: V ->           ; k73*V                                                     ;
    
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
    P    = 1.0     ;   
    RIGI = 1.0     ;
    IRF7 = 0.72205 ;
    V    = 6.9e-8  ;
'''

submodel_string = '''
    //Equations
    E2a: -> IFN         ; P*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)    ;
    E2b: IFN -> IFNe    ; k21*IFN                                           ;
    E4: -> STATP        ; P*k31*IFNe/(k32+k33*IFNe)                         ;
    E4b: STATP ->       ; t3*STATP                                          ;
    E5a: -> IRF7        ; P*k41*STATP+k42*IRF7P                             ;
    E5b: IRF7 ->        ; t4*IRF7                                           ; 
    E6a: -> IRF7P       ; P*k51*IRF7                                        ; 
    E6b: IRF7P ->       ; t5*IRF7P                                          ;
    E8a: -> V           ; P*(k71*V)/(1.0+k72*IFN*7E-5)                      ;
    E8b: V ->           ; k73*V                                             ;

    //Parameters
    k11 = 0.0       ;
    k12 = 9.746     ; 
    k13 = 12.511    ; 
    k14 = 13.562    ;
    k21 = 10.385    ;
    k31 = 45.922    ;
    k32 = 5.464     ;
    k33 = 0.068     ;
    t3  = 0.3       ;
    k41 = 0.115     ;
    k42 = 1.053     ;
    t4  = 0.3       ;
    k51 = 0.202     ;
    t5  = 0.3       ;
    k71 = 1.537     ;
    k72 = 47.883    ;
    k73 = 0.197     ;
    n   = 3.0       ;
    
    IFNe = 0.0      ;
    P = 1.0         ;

    //Initial Conditions
    IRF7 = 0.72205  ;
    V = 6.9e-8      ;
    RIGI = 1.0      ;
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

        # Changing initial values (if necessary)
        state = {}
        self.set_sbml_state(model_name='ODEModel', state=state)

    def step(self, mcs):
        #Step forward
        pass

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # set initial model parameters
        self.get_xml_element('IFNe_decay').cdata = self.sbml.ODEModel['t2'] * hours_to_mcs
        self.initial_infected = len(self.cell_list_by_type(self.I2))
        self.IFNe = self.sbml.ODEModel['IFNe']
        self.V = self.sbml.ODEModel['V']

        self.shared_steppable_vars['IFNe'] = self.V
        self.shared_steppable_vars['V'] = self.IFNe

        self.add_antimony_to_cell_types(model_string=submodel_string, model_name='submodel',cell_types=[self.I2],
                                        step_size=hours_to_mcs)

    def step(self, mcs):
        pIFNe = 0.0
        pV = 0.0
        secretor = self.get_field_secretor("IFNe")
        for cell in self.cell_list_by_type(self.I2):
            #Rule 3a
            k21 = self.sbml.ODEModel['k21'] * hours_to_mcs / self.initial_infected
            IFN = self.sbml.ODEModel['IFN']
            release = secretor.secreteInsideCellTotalCount(cell, k21 * IFN / cell.volume)
            pIFNe += abs(release.tot_amount)

            # Rule 7a
            k61 = self.sbml.ODEModel['k61'] * hours_to_mcs
            V = self.sbml.ODEModel['V']
            p_I2toDead = k61 * V
            if np.random.random() < p_I2toDead:
                cell.type = self.DEAD

        # Rule 3b
        t2 = self.sbml.ODEModel['t2'] * hours_to_mcs
        self.IFNe += pIFNe - t2 * self.IFNe
        self.shared_steppable_vars['IFNe'] = self.IFNe

        # Rule 8b
        k73 = self.sbml.ODEModel['k73'] * hours_to_mcs
        self.V += pV - k73 * self.V
        self.shared_steppable_vars['V'] = self.V

class PlotODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_infected = len(self.cell_list_by_type(self.I2))

        # Initialize Graphic Window for Amber Smith ODE model
        if (plot_ODEModel or plot_CellularizedModel):
            self.plot_win = self.add_new_plot_window(title='Jordan Model Cells',
                                                     x_axis_title='Hours',
                                                     y_axis_title='Number of Cells', x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False, config_options={'legend': True})

            self.plot_win2 = self.add_new_plot_window(title='Jordan Model IFNe',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Extracellular IFN', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            if plot_ODEModel:
                self.plot_win.add_plot("JP", style='Dots', color='red', size=5)
                self.plot_win2.add_plot("ODEVariable", style='Dots', color='blue', size=5)

            if plot_CellularizedModel:
                self.plot_win.add_plot("I2", style='Lines', color='red', size=5)
                self.plot_win2.add_plot("CC3DVariable", style='Lines', color='blue', size=5)

    def step(self, mcs):
        for cell in self.cell_list_by_type(self.I2):
            cell.sbml.submodel['IFNe'] = self.sbml.ODEModel['IFNe']
        self.timestep_sbml()

        self.V = 0.0
        for cell in self.cell_list_by_type(self.I2):
            self.V += cell.sbml.submodel['V']
        self.V /= self.initial_infected

        if plot_ODEModel:
            self.plot_win.add_data_point("JP", mcs * hours_to_mcs,self.sbml.ODEModel['P'])
            self.plot_win2.add_data_point("ODEVariable", mcs * hours_to_mcs, self.sbml.ODEModel['V'])

        if plot_CellularizedModel:
            num_I2 = len(self.cell_list_by_type(self.I2))
            self.plot_win.add_data_point("I2", mcs * hours_to_mcs, num_I2 / self.initial_infected)
            self.plot_win2.add_data_point("CC3DVariable", mcs * hours_to_mcs,self.V)