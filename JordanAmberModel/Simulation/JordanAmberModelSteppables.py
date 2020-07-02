from cc3d.core.PySteppables import *
import numpy as np
import os

plot_ODEModel = True
plot_CellModel = True
Data_writeout = False

## How to determine V
# 1 pulls from the scalar virus from the ODE original model (no feedback in the cellular model)
# 2 pulls from the scalar virus from the cellular model (feedback in the cellular model but no field)
# 3 pulls from the virus field
how_to_determine_V = 2

min_to_mcs = 10.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0 # hours/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 10.0  # 10 in the original model

'''Smith AP, Moquin DJ, Bernhauerova V, Smith AM. Influenza virus infection model with density dependence 
supports biphasic viral decay. Frontiers in microbiology. 2018 Jul 10;9:1554.'''

FluModel_string = '''        
        model FluModel()

        //State Variables and Transitions
        V1: -> T  ; -beta * V * T ;                             // Susceptible Cells
        V2: -> I1 ;  beta * V * T - k * I1 ;                    // Early Infected Cells
        V3: -> I2 ;  k * I1 - delta_d * I2 / (K_delta + I2) ;   // Late Infected Cells
        V4: -> V  ;  p * I2 - c * V ;                           // Extracellular Virus
        V5: -> D  ;  delta_d * I2 / (K_delta + I2) ;            // Cleared Infected Cells (for Bookkeeping)

        //Parameters
        beta = 2.4* 10^(-4) ;                                   // Virus Infective
        p = 1.6 ;                                               // Virus Production
        c = 13.0 ;                                              // Virus Clearance
        k = 4.0 ;                                               // Eclipse phase
        delta_d = 1.6 * 10^6 ;                                  // Infected Cell Clearance
        K_delta = 4.5 * 10^5 ;                                  // Half Saturation Constant         

        // Initial Conditions ;
        T0 = 1.0*10^7;
        T = T0  ;                                               // Initial Number of Uninfected Cells
        I1 = 75.0 ;                                             // Initial Number of Infected Cells
end'''

'''Jordan J. A. Weaver and Jason E. Shoemaker. Mathematical Modeling of RNA Virus Sensing Pathways Reveal Paracrine Signaling as the Primary Factor 
Regulating Excessive Cytokine Production'''

# Modified code to make IFNe production dependent on the number of infected cells (P)
IFNModel_string = '''
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
    t4  = 0.75      ;
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
    IRF7 = 0.0      ;
    V    = 6.9e-8   ;
'''

# Viral Replication Model
viral_model_string = '''
    E7a: P ->           ; P*k61*V                    ;
    E8a: -> V           ; P*k71*V/(1.0+k72*IFN*7E-5) ;
    E8b: V ->           ; k73*V                      ;

    //Parameters
    k61 = 0.635     ;
    k71 = 1.537     ;
    k72 = 47.883    ;
    k73 = 0.197     ;

    //Initial Conditions
    P    =  1.0     ;
    V    =  0.0     ; 
    
    //Inputs
    IFN  =  0.0     ;
'''

# Modified code to make IFNe production dependent on the number of infected cells (P)
IFN_model_string = '''
    //Equations
    E2a: -> IFN         ; P*(k12*(V^n)/(k13+(V^n))+k14*IRF7P)               ;
    E2b: IFN ->         ; k21*IFN                                           ;
    E4a: -> STATP       ; P*k31*IFNe/(k32+k33*IFNe)                         ;
    E4b: STATP ->       ; t3*STATP                                          ;
    E5a: -> IRF7        ; P*(k41*STATP+k42*IRF7P)                           ;
    E5b: IRF7 ->        ; t4*IRF7                                           ;
    E6a: -> IRF7P       ; P*k51*IRF7                                        ;
    E6b: IRF7P ->       ; t5*IRF7P                                          ;

    //Parameters
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
    n   = 3.0       ;

    // Inputs
    P    = 0.0      ;
    V    = 0.0      ;
    IFNe = 0.0      ;
'''

class ODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Uptading max simulation steps using scaling factor to simulate 10 days
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs

        # Adding free floating antimony model
        self.add_free_floating_antimony(model_string=FluModel_string, model_name='FluModel',
                                        step_size=days_to_mcs)

        self.add_free_floating_antimony(model_string=IFNModel_string, model_name='IFNModel',
                                        step_size=hours_to_mcs)

        self.add_antimony_to_cell_types(model_string=viral_model_string, model_name='VModel',
                                        cell_types=[self.U], step_size=hours_to_mcs)

        self.add_antimony_to_cell_types(model_string=IFN_model_string, model_name='IModel',
                                        cell_types=[self.U], step_size=hours_to_mcs)

        # Changing initial values according to discussions with Amber Smith
        self.sbml.FluModel['I1'] = 0.0
        self.sbml.FluModel['V'] = 75.0

    def step(self, mcs):
        self.timestep_sbml()

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # set initial model parameters
        self.initial_uninfected = len(self.cell_list_by_type(self.U))
        self.ExtracellularVirus = self.sbml.FluModel['V']
        self.get_xml_element('virus_decay').cdata = self.sbml.FluModel['c'] * days_to_mcs
        self.ExtracellularIFN = self.sbml.IFNModel['IFNe']
        self.get_xml_element('IFNe_decay').cdata = self.sbml.IFNModel['k73'] * hours_to_mcs

    def step(self, mcs):
        secretorV = self.get_field_secretor("Virus")
        secretorIFN = self.get_field_secretor("IFNe")

        for cell in self.cell_list_by_type(self.U):
            ## U to I1 transition - Amber Model
            # Determine V from scalar virus from the ODE
            if how_to_determine_V == 1:
                b = self.sbml.FluModel['beta'] * self.sbml.FluModel['T0'] * days_to_mcs
                V = self.sbml.FluModel['V'] / self.sbml.FluModel['T0']
            # Determine V from scalar virus from the cellular model
            if how_to_determine_V == 2:
                b = self.sbml.FluModel['beta'] * self.initial_uninfected * days_to_mcs
                V = self.ExtracellularVirus / self.initial_uninfected
            # Determine V from the virus field
            if how_to_determine_V == 3:
                b = self.sbml.FluModel['beta'] * self.initial_uninfected * days_to_mcs
                V = secretorV.amountSeenByCell(cell)
            # V1: T -> U ; beta * V * T
            p_UtoI1 = b * V
            if np.random.random() < p_UtoI1:
                cell.type = self.I1
                cell.sbml.VModel['V'] = 6.9e-8

        ## I1 to I2 transition - Amber Model
        # V2: I1 -> I2 ;k * I1
        k = self.sbml.FluModel['k'] * days_to_mcs
        p_T1oI2 = k
        for cell in self.cell_list_by_type(self.I1):
            if np.random.random() < p_T1oI2:
                cell.type = self.I2

        ## P to D transition - Jordan Model
        # E7a: P -> ; P * k61 * V;
        for cell in self.cell_list_by_type(self.I1, self.I2):
            k61 = cell.sbml.VModel['k61'] * hours_to_mcs
            V = cell.sbml.VModel['V']
            P = cell.sbml.VModel['P']
            p_I2toD = k61 * V * (1-P)
            if np.random.random() < p_I2toD:
                cell.type = self.DEAD

        ## Updating values of intracellular models
        for cell in self.cell_list_by_type(self.U,self.I1, self.I2):
            ## Inputs to the INF model
            cell.sbml.IModel['V'] = cell.sbml.VModel['V']
            cell.sbml.IModel['P'] = cell.sbml.VModel['P']
            IFNe = secretorIFN.amountSeenByCell(cell)
            cell.sbml.IModel['IFNe'] = IFNe # Rescale
            ## Inputs to the Virus model
            cell.sbml.VModel['IFN'] = cell.sbml.IModel['IFN']

        ## Production of extracellular virus - Jordan Model
        # E8b: V -> ; k73 * V
        V = self.ExtracellularVirus
        k73 = self.sbml.IFNModel['k73'] * hours_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            Virus = cell.sbml.VModel['V']
            p = k73 * Virus * self.sbml.FluModel['T0'] #Rescale?
            release = secretorV.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.ExtracellularVirus += release.tot_amount
        c = self.sbml.FluModel['c'] * days_to_mcs
        self.ExtracellularVirus -= c * V

        ## Measure amount of extracellular virus field
        self.ExtracellularVirus_Field = 0
        for cell in self.cell_list:
            V = secretorV.amountSeenByCell(cell)
            self.ExtracellularVirus_Field += V

        ## Production of extracellular IFN - Jordan Model
        # E2b: IFN -> IFNe; k21 * IFN ;
        I = self.ExtracellularIFN
        k21 = self.sbml.IFNModel['k21'] * hours_to_mcs
        for cell in self.cell_list_by_type(self.I1,self.I2):
            intracellularIFN = cell.sbml.IModel['IFN']
            p = k21 * intracellularIFN #/ 3.0 #rescale?
            release = secretorIFN.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.ExtracellularIFN += release.tot_amount
        # E3a: IFNe -> ; t2*IFNe ;
        t2 = self.sbml.IFNModel['t2'] * hours_to_mcs
        self.ExtracellularIFN -= t2 * I

        ## Measure amount of extracellular IFN field
        self.ExtracellularIFN_Field = 0
        for cell in self.cell_list:
            I = secretorIFN.amountSeenByCell(cell)
            self.ExtracellularIFN_Field += I

        # Dictonary to pass information between steppables
        self.shared_steppable_vars['ExtracellularVirus'] = self.ExtracellularVirus
        self.shared_steppable_vars['ExtracellularVirus_Field'] = self.ExtracellularVirus_Field
        self.shared_steppable_vars['ExtracellularIFN'] = self.ExtracellularIFN
        self.shared_steppable_vars['ExtracellularIFN_Field'] = self.ExtracellularIFN_Field

class FluPlotSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_uninfected = len(self.cell_list_by_type(self.U))
        if (plot_ODEModel == True) or (plot_CellModel == True):
            self.plot_win = self.add_new_plot_window(title='Flu Model Cells',
                                                     x_axis_title='Hours',
                                                     y_axis_title='Variables', x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False, config_options={'legend': True})

            self.plot_win2 = self.add_new_plot_window(title='Flu Model Virus',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Virus', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            if plot_ODEModel == True:
                self.plot_win.add_plot("ODET", style='Dots', color='blue', size=5)
                self.plot_win.add_plot("ODEI1", style='Dots', color='orange', size=5)
                self.plot_win.add_plot("ODEI2", style='Dots', color='red', size=5)
                self.plot_win.add_plot("ODED", style='Dots', color='purple', size=5)
                self.plot_win2.add_plot("ODEV", style='Dots', color='blue', size=5)

            if plot_CellModel == True:
                self.plot_win.add_plot("CC3DT", style='Lines', color='blue', size=5)
                self.plot_win.add_plot("CC3DI1", style='Lines', color='orange', size=5)
                self.plot_win.add_plot("CC3DI2", style='Lines', color='red', size=5)
                self.plot_win.add_plot("CC3DD", style='Lines', color='purple', size=5)
                self.plot_win2.add_plot("CC3DV", style='Lines', color='blue', size=5)

    def step(self, mcs):
        if (plot_ODEModel == True) or (plot_CellModel == True):
            if plot_ODEModel == True:
                self.plot_win.add_data_point("ODET", mcs * days_to_mcs * 24.0,
                                             self.sbml.FluModel['T'] / self.sbml.FluModel['T0'])
                self.plot_win.add_data_point("ODEI1", mcs * days_to_mcs * 24.0,
                                             self.sbml.FluModel['I1'] / self.sbml.FluModel['T0'])
                self.plot_win.add_data_point("ODEI2", mcs * days_to_mcs * 24.0,
                                             self.sbml.FluModel['I2'] / self.sbml.FluModel['T0'])
                self.plot_win.add_data_point("ODED", mcs * days_to_mcs * 24.0,
                                             self.sbml.FluModel['D'] / self.sbml.FluModel['T0'])
                self.plot_win2.add_data_point("ODEV", mcs * days_to_mcs * 24.0,
                                              np.log10(self.sbml.FluModel['V']))

            if plot_CellModel == True:
                self.plot_win.add_data_point("CC3DT", mcs * days_to_mcs * 24.0,
                                             len(self.cell_list_by_type(self.U)) / self.initial_uninfected)
                self.plot_win.add_data_point("CC3DI1", mcs * days_to_mcs * 24.0,
                                             len(self.cell_list_by_type(self.I1)) / self.initial_uninfected)
                self.plot_win.add_data_point("CC3DI2", mcs * days_to_mcs * 24.0,
                                             len(self.cell_list_by_type(self.I2)) / self.initial_uninfected)
                self.plot_win.add_data_point("CC3DD", mcs * days_to_mcs * 24.0,
                                             len(self.cell_list_by_type(self.DEAD)) / self.initial_uninfected)
                self.plot_win2.add_data_point("CC3DV", mcs * days_to_mcs * 24.0,
                                              np.log10(self.shared_steppable_vars['ExtracellularVirus_Field']))

class IFNPlotSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_uninfected = len(self.cell_list_by_type(self.U))
        # Initialize Graphic Window for Jordan IFN model
        if (plot_ODEModel == True) or (plot_CellModel == True):
            self.plot_win3 = self.add_new_plot_window(title='P',
                                                     x_axis_title='Hours',
                                                     y_axis_title='Number of Cells', x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False, config_options={'legend': True})

            self.plot_win5 = self.add_new_plot_window(title='IFNe',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Extracellular IFN', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win8 = self.add_new_plot_window(title='V',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            if plot_ODEModel:
                self.plot_win3.add_plot("ODEP", style='Dots', color='red', size=5)
                self.plot_win3.add_plot("FLUP", style='Dots', color='orange', size=5)
                self.plot_win5.add_plot("ODEIFNe", style='Dots', color='purple', size=5)
                self.plot_win8.add_plot("ODEV", style='Dots', color='yellow', size=5)

            if plot_CellModel:
                self.plot_win3.add_plot("CC3DP", style='Lines', color='red', size=5)
                self.plot_win5.add_plot("CC3DIFNe", style='Lines', color='purple', size=5)
                self.plot_win8.add_plot("CC3DV", style='Lines', color='yellow', size=5)

    def step(self, mcs):
        if plot_ODEModel:
            FLUP = (self.sbml.FluModel['T'] + self.sbml.FluModel['I1'] + self.sbml.FluModel['I2'] ) / self.sbml.FluModel['T0']
            self.plot_win3.add_data_point("ODEP", mcs * hours_to_mcs,self.sbml.IFNModel['P'])
            self.plot_win3.add_data_point("FLUP", mcs * hours_to_mcs, FLUP)
            self.plot_win5.add_data_point("ODEIFNe", mcs * hours_to_mcs, self.sbml.IFNModel['IFNe'])
            self.plot_win8.add_data_point("ODEV", mcs * hours_to_mcs, self.sbml.IFNModel['V'])

        if plot_CellModel:
            P = len(self.cell_list_by_type(self.U,self.I1,self.I2))/self.initial_uninfected
            maxV = 0.0
            L = len(self.cell_list_by_type(self.I1,self.I2))
            for cell in self.cell_list_by_type(self.I1,self.I2):
                maxV += cell.sbml.VModel['V'] / L
            self.plot_win3.add_data_point("CC3DP", mcs * hours_to_mcs, P)
            self.plot_win5.add_data_point("CC3DIFNe", mcs * hours_to_mcs, self.shared_steppable_vars['ExtracellularIFN_Field']/self.initial_uninfected)
            self.plot_win8.add_data_point("CC3DV", mcs * hours_to_mcs, maxV)

class Data_OutputSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        if Data_writeout:
            folder_path = '/Users/Josua/AmberFluModelData/'
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)

            file_name = 'ODEModel.txt'
            self.output1 = open(folder_path + file_name, 'w')
            self.output1.write(
                "%s,%s,%s,%s,%s,%s\n" % ('Time', 'T', 'I1', 'I2', 'D', 'V'))
            self.output1.flush()

            file_name = 'CC3DModel.txt'
            self.output2 = open(folder_path + file_name, 'w')
            self.output2.write(
                "%s,%s,%s,%s,%s,%s\n" % ('Time', 'T', 'I1', 'I2', 'D', 'V'))
            self.output2.flush()

    def step(self, mcs):
        if Data_writeout:
            # Record variables from ODE model
            d = mcs * days_to_mcs
            AT = self.sbml.FluModel['T']
            AI1 = self.sbml.FluModel['I1']
            AI2 = self.sbml.FluModel['I2']
            AD = self.sbml.FluModel['D']
            AV = self.sbml.FluModel['V']

            self.output1.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (d, AT, AI1, AI2, AD, AV))
            self.output1.flush()

            # Record variables from Cellularized Model
            U = len(self.cell_list_by_type(self.U))
            I1 = len(self.cell_list_by_type(self.I1))
            I2 = len(self.cell_list_by_type(self.I2))
            D = len(self.cell_list_by_type(self.DEAD))
            V = self.shared_steppable_vars['ExtracellularVirus_Field']

            self.output2.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (d, U, I1, I2, D, V))
            self.output2.flush()

    def finish(self):
        if Data_writeout:
            self.output1.close()
            self.output2.close()