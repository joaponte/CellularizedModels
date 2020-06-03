from cc3d.core.PySteppables import *
import numpy as np
import os

min_to_mcs = 5.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0  # hrs/mcs
hours_to_simulate = 36.0

plot_ODEModel = True
plot_CellularizedModel =  True

'''Jordan J. A. Weaver and Jason E. Shoemaker. Mathematical Modeling of RNA Virus Sensing Pathways Reveal Paracrine Signaling as the Primary Factor 
Regulating Excessive Cytokine Production'''

model_string = '''
    //Equations
    E2: -> IFN   ; P*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)-k21*IFN   ; // Intracellular IFN
    E3: -> IFNe  ; k21*IFN-t2*IFNe                                        ; // Extracellular IFN
    E4: -> STATP ; P*k31*IFNe/(k32+k33*IFNe)-t3*STATP                       ; // Intracellular STATP
    E5: -> IRF7  ; P*(k41*STATP+k42*IRF7P)-t4*IRF7                          ; // Intracellular IRF7
    E6: -> IRF7P ; P*k51*IRF7-t5*IRF7P                                      ; // Intracellular IRF7P
    E7: -> P     ; - k61*P*V                                                ; // Infected Cells
    E8: -> V     ; P*k71*V/(1+k72*IFNe)-k73*V                               ; // Intracellular Virus

    //Parameters
    k11 = 0.0       ;  // 
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
    k72 = 0.0       ; //47.883
    k73 = 0.197     ;
    n   = 3.0       ;

    //Initial Conditions
    P    = 1.0     ;   
    RIGI = 1.0     ;
    IRF7 = 0.72205 ;
    V    = 6.9e-8  ;
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
        self.timestep_sbml()

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # set initial model parameters
        self.initial_infected = len(self.cell_list_by_type(self.I2))
        self.IFNe = self.sbml.ODEModel['IFNe']
        self.V = self.sbml.ODEModel['V']
        self.get_xml_element('IFNe_decay').cdata = self.sbml.ODEModel['t2'] * hours_to_mcs

        for cell in self.cell_list_by_type(self.I2):
            cell.dict['IFN'] = self.sbml.ODEModel['IFN']
            cell.dict['STATP'] = self.sbml.ODEModel['STATP']
            cell.dict['IRF7'] = self.sbml.ODEModel['IRF7']
            cell.dict['IRF7P'] = self.sbml.ODEModel['IRF7P']
            cell.dict['V'] = self.sbml.ODEModel['V']

    def step(self, mcs):
        secretor = self.get_field_secretor("IFNe")
        for cell in self.cell_list_by_type(self.I2):
            # Rule E6
            p_I2toDead = self.sbml.ODEModel['k61'] * self.sbml.ODEModel['V'] * hours_to_mcs
            if np.random.random() < p_I2toDead:
                cell.type = self.DEAD

            #Rule E3
            k21 = self.sbml.ODEModel['k21'] * hours_to_mcs
            IFN = self.sbml.ODEModel['IFN'] / self.initial_infected
            release = secretor.secreteInsideCellTotalCount(cell, k21 * IFN / cell.volume)
            self.IFNe += release.tot_amount
        self.IFNe -= self.sbml.ODEModel['t2'] * self.IFNe * hours_to_mcs

class PlotODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_infected = len(self.cell_list_by_type(self.I2))
        self.IFNe = self.sbml.ODEModel['IFNe']

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
                self.plot_win2.add_plot("JIFNe", style='Dots', color='blue', size=5)

            if plot_CellularizedModel:
                self.plot_win.add_plot("I2", style='Lines', color='red', size=5)
                self.plot_win2.add_plot("IFNe", style='Lines', color='blue', size=5)

    def step(self, mcs):
        k21 = self.sbml.ODEModel['k21'] * hours_to_mcs
        IFN = self.sbml.ODEModel['IFN'] / self.initial_infected
        for cell in self.cell_list_by_type(self.I2):
            self.IFNe += k21 * IFN
        self.IFNe -= self.sbml.ODEModel['t2'] * self.IFNe * hours_to_mcs

        if plot_ODEModel:
            self.plot_win.add_data_point("JP", mcs * hours_to_mcs,
                                     self.sbml.ODEModel['P'])
            self.plot_win2.add_data_point("JIFNe", mcs * hours_to_mcs,
                                     self.sbml.ODEModel['IFNe'])

        if plot_CellularizedModel:
            num_I2 = len(self.cell_list_by_type(self.I2))
            self.plot_win.add_data_point("I2", mcs * hours_to_mcs, num_I2/self.initial_infected)
            self.plot_win2.add_data_point("IFNe", mcs * hours_to_mcs, self.IFNe)