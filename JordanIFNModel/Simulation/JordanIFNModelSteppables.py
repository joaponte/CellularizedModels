from cc3d.core.PySteppables import *
import numpy as np
import os

min_to_mcs = 1.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0  # hrs/mcs
hours_to_simulate = 36.0

'''Jordan J. A. Weaver and Jason E. Shoemaker. Mathematical Modeling of RNA Virus Sensing Pathways Reveal Paracrine Signaling as the Primary Factor 
Regulating Excessive Cytokine Production'''

model_string = '''
    //Equations
    E2: -> IFN   ; P*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)-k21*IFN   ; // Intracellular IFN
    E3: -> IFNe  ; k21*IFN-t2*IFNe                                          ; // Extracellular IFN
    E4: -> STATP ; P*k31*IFNe/(k32+k33*IFNe)-t3*STATP                       ; // Intracellular STATP
    E5: -> IRF7  ; P*(k41*STATP+k42*IRF7P)-t4*IRF7                          ; // Intracellular IRF7
    E6: -> IRF7P ; P*k51*IRF7-t5*IRF7P                                      ; // Intracellular IRF7P
    E7: -> P     ; - k61*P*V                                                ; // Infected Cells
    E8: -> V     ; k71*P*V/(1+k72*IFNe)-k73*V                               ; // Intracellular Virus

    //Parameters
    k11 = 0.0    ;  // 
    k12 = 9.746  ; 
    k13 = 12.511 ; 
    k14 = 13.562 ;
    k21 = 10.385 ;
    t2  = 3.481  ;
    k31 = 45.922 ;
    k32 = 5.464  ;
    k33 = 0.068  ;
    t3  = 0.3    ;
    k41 = 0.115  ;
    k42 = 1.053  ;
    t4  = 0.3    ;
    k51 = 0.202  ;
    t5  = 0.3    ;
    k61 = 0.635  ;
    k71 = 1.537  ;
    k72 = 47.883 ;
    k73 = 0.197  ;
    n   = 3.0    ;

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
        self.FNe = self.sbml.ODEModel['IFNe']
        self.get_xml_element('IFNe_decay').cdata = self.sbml.JordansModel['t2'] * days_to_mcs * 24.0

        for cell in self.cell_list_by_type(self.I2):
            cell.dict['IFN'] = self.sbml.ODEModel['IFN']
            cell.dict['STATP'] = self.sbml.ODEModel['STATP']
            cell.dict['IRF7'] = self.sbml.ODEModel['IRF7']
            cell.dict['IRF7P'] = self.sbml.ODEModel['IRF7P']
            cell.dict['V'] = self.sbml.ODEModel['V']

    def step(self, mcs):
        # Transtion from Infection to Dead
        kd = 1.0 * mcs
        for cell in self.cell_list_by_type(self.Infected):
            if cell.dict['Time_since_infection'] > 18.0:
                cell.type = self.DEAD

        # Transition rule to infected cell
        for cell in self.cell_list_by_type(self.Infected):
                cell.dict['Time_since_infection'] = mcs

        p_Infected_to_Dead = self.sbml.JordansModel['k61'] /self.initial_uninfected * self.sbml.JordansModel['V'] * days_to_mcs
        for cell in self.cell_list_by_type(self.Infected):
            if np.random.random() < p_UtoI1:
                cell.type = self.DEAD

        # Transition rule from I1 to I2
        k = self.sbml.ambersmithsimple['k'] * days_to_mcs
        p_T1oI2 = k
        for cell in self.cell_list_by_type(self.I1):
            if np.random.random() < p_T1oI2:
                cell.type = self.I2

        # Transition rule from I2 to D
        K_delta = self.sbml.ambersmithsimple['K_delta'] / self.sbml.ambersmithsimple['T0'] * self.initial_uninfected
        delta_d = self.sbml.ambersmithsimple['delta_d'] / self.sbml.ambersmithsimple['T0'] * self.initial_uninfected
        I2 = len(self.cell_list_by_type(self.I2))
        p_T2toD = delta_d / (K_delta + I2) * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            if np.random.random() < p_T2toD:
                cell.type = self.DEAD

        # Production of extracellular virus
        secretor = self.get_field_secretor("Virus")
        V = self.ExtracellularVirus
        p = self.sbml.ambersmithsimple['p'] / self.initial_uninfected * self.sbml.ambersmithsimple['T0'] * days_to_mcs
        c = self.sbml.ambersmithsimple['c'] * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            release = secretor.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.ExtracellularVirus += release.tot_amount
        self.ExtracellularVirus -= c * V

        # Measure amount of extracellular virus field
        self.ExtracellularVirus_Field = 0
        for cell in self.cell_list:
            uptake_probability = 0.0000001
            uptake = secretor.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            V = abs(uptake.tot_amount) / uptake_probability
            self.ExtracellularVirus_Field += V
            secretor.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)