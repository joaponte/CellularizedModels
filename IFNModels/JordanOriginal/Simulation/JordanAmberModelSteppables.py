from cc3d.core.PySteppables import *
import numpy as np
import os

plot_ODEModel = True
plot_CellModel = True
OutputData = True

how_to_determine_IFNe = 2 # Determines the IFNe from the ODE model (1) from Cell model as scalar (2) or from field (3)

min_to_mcs = 10.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0 # hours/mcs
hours_to_simulate = 30.0

IFNe_diffusion_coefficient = 1.0/10.0 #vl^2 / min

'''Jordan J. A. Weaver and Jason E. Shoemaker. Mathematical Modeling of RNA Virus Sensing Pathways Reveal Paracrine Signaling as the Primary Factor 
Regulating Excessive Cytokine Production'''

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
    E8a: -> V           ; P*(k71*V)/(1.0+k72*IFNe*7E-5)                     ;
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
    E7a: H ->           ; H*k61*V                     ;
    E8a: -> V           ; H*k71*V/(1.0+k72*IFNe*7E-5) ;
    E8b: V ->           ; k73*V                       ;

    //Parameters
    k61 = 0.635     ;
    k71 = 1.537     ;
    k72 = 47.883    ;
    k73 = 0.197     ;

    //Initial Conditions
    V =  6.9e-8      ; 
    H = 1.0          ;
    
    //Inputs
    IFNe  =  0.0     ;
'''

IFN_model_string = '''
    //Equations
    E2a: -> IFN         ; H*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)    ;
    E2b: IFN ->         ; k21*IFN                                           ;
    E4a: -> STATP       ; H*k31*IFNe/(k32+k33*IFNe)                         ;
    E4b: STATP ->       ; t3*STATP                                          ;
    E5a: -> IRF7        ; H*(k41*STATP+k42*IRF7P)                           ;
    E5b: IRF7 ->        ; t4*IRF7                                           ;
    E6a: -> IRF7P       ; H*k51*IRF7                                        ;
    E6b: IRF7P ->       ; t5*IRF7P                                          ;

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
    t4  = 0.75      ;
    k51 = 0.202     ;
    t5  = 0.3       ;
    n   = 3.0       ;
    RIGI = 1.0      ;
    
    // Inputs
    H    = 0.0      ;
    IFNe = 0.0      ;
    V = 0.0         ;
'''

class ODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Store Initial Number of Cells
        self.shared_steppable_vars['InitialNumberCells'] = len(self.cell_list)

        # Set Max Simulation Steps
        self.get_xml_element('simulation_steps').cdata = hours_to_simulate / hours_to_mcs

        # Load Original ODE Model
        self.add_free_floating_antimony(model_string=IFNModel_string, model_name='IFNModel',
                                        step_size=hours_to_mcs)

        # Load Viral Model inside Cells
        self.add_antimony_to_cell_types(model_string=viral_model_string, model_name='VModel',
                                        cell_types=[self.I2], step_size=hours_to_mcs)

        # Load IFN Model inside Cells
        self.add_antimony_to_cell_types(model_string=IFN_model_string, model_name='IModel',
                                        cell_types=[self.I2], step_size=hours_to_mcs)

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Set IFNe diffusion parameters
        self.get_xml_element('IFNe_dc').cdata = IFNe_diffusion_coefficient * min_to_mcs
        self.get_xml_element('IFNe_decay').cdata = self.sbml.IFNModel['t2'] * hours_to_mcs
        self.shared_steppable_vars['ExtracellularIFN_Scalar'] = self.sbml.IFNModel['IFNe']

        # Set secretors
        self.secretorIFN = self.get_field_secretor("IFNe")

        # Assign cell lifetime
        for cell in self.cell_list:
            cell.dict['lifetime'] = np.random.normal(24,2)

    def step(self, mcs):
        ## Measure amount of IFNe in the Field
        self.shared_steppable_vars['ExtracellularIFN_Field'] = 0
        for cell in self.cell_list_by_type(self.I2):
            self.shared_steppable_vars['ExtracellularIFN_Field'] += self.secretorIFN.amountSeenByCell(cell)

        ## Production of IFNe
        # E2b: IFN -> IFNe; k21 * IFN ;
        self.total_IFNe_production = 0.0
        k21 = self.sbml.IFNModel['k21'] * hours_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            intracellularIFN = cell.sbml.IModel['IFN']
            p = k21 * intracellularIFN
            release = self.secretorIFN.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.total_IFNe_production += release.tot_amount

        ## Decay of IFNe
        # E3a: IFNe -> ; t2*IFNe ;
        I = self.shared_steppable_vars['ExtracellularIFN_Scalar']
        t2 = self.sbml.IFNModel['t2'] * hours_to_mcs
        self.total_decay = t2 * I
        ## Update Scalar IFNe
        self.shared_steppable_vars['ExtracellularIFN_Scalar'] += self.total_IFNe_production - self.total_decay

        ## P to D transition
        # E7a: P -> ; P * k61 * V;
        for cell in self.cell_list_by_type(self.I2):
            k61 = cell.sbml.VModel['k61'] * hours_to_mcs
            H = cell.sbml.VModel['H']
            V = cell.sbml.VModel['V']
            r = k61 * V * (1-H)
            p_I2toD = 1.0 - np.exp(-r)
            if np.random.random() < p_I2toD:
                cell.type = self.DEAD

        ## Addtiional P to D transition
        # E7a: P -> ; P * k61 * V;
        for cell in self.cell_list:
            cell.dict['lifetime'] -= hours_to_mcs
            if cell.dict['lifetime'] <= 0.0:
                cell.type = self.DEAD

        ## Updating Cellular Models
        L = len(self.cell_list_by_type(self.I2))
        for cell in self.cell_list:
            if how_to_determine_IFNe == 1:
                cell.sbml.VModel['IFNe'] = self.sbml.IFNModel['IFNe']
                cell.sbml.IModel['IFNe'] = self.sbml.IFNModel['IFNe']
            if how_to_determine_IFNe == 2:
                cell.sbml.VModel['IFNe'] = self.shared_steppable_vars['ExtracellularIFN_Scalar'] / self.shared_steppable_vars['InitialNumberCells']
                cell.sbml.IModel['IFNe'] = self.shared_steppable_vars['ExtracellularIFN_Scalar'] / self.shared_steppable_vars['InitialNumberCells']
            if how_to_determine_IFNe == 3:
                cell.sbml.VModel['IFNe'] = self.secretorIFN.amountSeenByCell(cell)
                cell.sbml.IModel['IFNe'] = self.secretorIFN.amountSeenByCell(cell)
            cell.sbml.IModel['H'] = cell.sbml.VModel['H']
            cell.sbml.IModel['V'] = cell.sbml.VModel['V']

        self.timestep_sbml()

class IFNPlotSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_infected = len(self.cell_list_by_type(self.U))
        # Initialize Graphic Window for Jordan IFN model
        if (plot_ODEModel == True) or (plot_CellModel == True):
            self.plot_win1 = self.add_new_plot_window(title='V',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win2 = self.add_new_plot_window(title='H',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})


            self.plot_win3 = self.add_new_plot_window(title='P',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win4 = self.add_new_plot_window(title='IFNe',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win5 = self.add_new_plot_window(title='STATP',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win6 = self.add_new_plot_window(title='IRF7',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win7 = self.add_new_plot_window(title='IRF7P',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            self.plot_win8 = self.add_new_plot_window(title='INF',
                                                      x_axis_title='Hours',
                                                      y_axis_title='Variable', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False, config_options={'legend': True})

            if plot_ODEModel:
                self.plot_win1.add_plot("ODEV", style='Dots', color='yellow', size=5)
                self.plot_win2.add_plot("ODEH", style='Dots', color='white', size=5)
                self.plot_win3.add_plot("ODEP", style='Dots', color='red', size=5)
                self.plot_win4.add_plot("ODEIFNe", style='Dots', color='orange', size=5)
                self.plot_win5.add_plot("ODESTATP", style='Dots', color='blue', size=5)
                self.plot_win6.add_plot("ODEIRF7", style='Dots', color='green', size=5)
                self.plot_win7.add_plot("ODEIRF7P", style='Dots', color='purple', size=5)
                self.plot_win8.add_plot("ODEIFN", style='Dots', color='magenta', size=5)

            if plot_CellModel:
                self.plot_win1.add_plot("CC3DV", style='Lines', color='yellow', size=5)
                self.plot_win2.add_plot("CC3DH", style='Lines', color='white', size=5)
                self.plot_win3.add_plot("CC3DP", style='Lines', color='red', size=5)
                self.plot_win4.add_plot("CC3DIFNe_Scalar", style='Lines', color='orange', size=5)
                self.plot_win4.add_plot("CC3DIFNe_Field", style='Lines', color='yellow', size=5)
                self.plot_win5.add_plot("CC3DSTATP", style='Lines', color='blue', size=5)
                self.plot_win6.add_plot("CC3DIRF7", style='Lines', color='green', size=5)
                self.plot_win7.add_plot("CC3DIRF7P", style='Lines', color='purple', size=5)
                self.plot_win8.add_plot("CC3DIFN", style='Lines', color='magenta', size=5)

    def step(self, mcs):
        if plot_ODEModel:
            P = len(self.cell_list_by_type(self.I2))/self.shared_steppable_vars['InitialNumberCells']
            self.plot_win1.add_data_point("ODEV", mcs * hours_to_mcs, self.sbml.IFNModel['V'])
            self.plot_win2.add_data_point("ODEH", mcs * hours_to_mcs, self.sbml.IFNModel['P'])
            self.plot_win3.add_data_point("ODEP", mcs * hours_to_mcs, self.sbml.IFNModel['P'])
            self.plot_win4.add_data_point("ODEIFNe", mcs * hours_to_mcs, self.sbml.IFNModel['IFNe'] * P)
            self.plot_win5.add_data_point("ODESTATP", mcs * hours_to_mcs, self.sbml.IFNModel['STATP'])
            self.plot_win6.add_data_point("ODEIRF7", mcs * hours_to_mcs, self.sbml.IFNModel['IRF7'])
            self.plot_win7.add_data_point("ODEIRF7P", mcs * hours_to_mcs, self.sbml.IFNModel['IRF7P'])
            self.plot_win8.add_data_point("ODEIFN", mcs * hours_to_mcs, self.sbml.IFNModel['IFN'])

        if plot_CellModel:
            L = len(self.cell_list_by_type(self.I2))
            avgV = 0.0
            avgH = 0.0
            avgSTATP = 0.0
            avgIRF7 = 0.0
            avgIRF7P = 0.0
            avgIFN = 0.0
            for cell in self.cell_list_by_type(self.I2):
                avgV += cell.sbml.VModel['V'] / L
                avgH += cell.sbml.VModel['H'] / L
                avgSTATP += cell.sbml.IModel['STATP'] / L
                avgIRF7 += cell.sbml.IModel['IRF7'] / L
                avgIRF7P += cell.sbml.IModel['IRF7P'] / L
                avgIFN += cell.sbml.IModel['IFN'] / L

            self.plot_win1.add_data_point("CC3DV", mcs * hours_to_mcs, avgV)
            self.plot_win2.add_data_point("CC3DH", mcs * hours_to_mcs, avgH)
            self.plot_win3.add_data_point("CC3DP", mcs * hours_to_mcs,
                                          L / self.shared_steppable_vars['InitialNumberCells'])
            self.plot_win4.add_data_point("CC3DIFNe_Scalar", mcs * hours_to_mcs,
                                          self.shared_steppable_vars['ExtracellularIFN_Scalar']
                                          / self.shared_steppable_vars['InitialNumberCells'])
            self.plot_win4.add_data_point("CC3DIFNe_Field", mcs * hours_to_mcs,
                                          self.shared_steppable_vars['ExtracellularIFN_Field']
                                          / self.shared_steppable_vars['InitialNumberCells'])
            self.plot_win5.add_data_point("CC3DSTATP", mcs * hours_to_mcs, avgSTATP)
            self.plot_win6.add_data_point("CC3DIRF7", mcs * hours_to_mcs, avgIRF7)
            self.plot_win7.add_data_point("CC3DIRF7P", mcs * hours_to_mcs, avgIRF7P)
            self.plot_win8.add_data_point("CC3DIFN", mcs * hours_to_mcs, avgIFN)

class OutputSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        if OutputData:
            folder_path = '/Users/Josua/Data/'
            # folder_path = '/N/u/joaponte/Carbonate/FluModel/Output/'
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            Replicate = 1
            # Output ODE Data
            file_name1 = 'JordanOriginalODE_%i.txt' % Replicate
            self.output1 = open(folder_path + file_name1, 'w')
            self.output1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %
                              ('Time','ODEV', 'ODEH', 'ODEP', 'ODEIFNe','ODESTATP','ODEIRF7','ODEIRF7P','ODEIFN'))
            self.output1.flush()

            # Output CC3D Data
            file_name2 = 'JordanOriginalCC3D_%i.txt' % Replicate
            self.output2 = open(folder_path + file_name2, 'w')
            self.output2.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %
                              ('Time','CC3DV','CC3DH','CC3DP','CC3DIFNe_Scalar','CC3DIFNe_Field','CC3DSTATP','CC3DIRF7',
                               'CC3DIRF7P','CC3DIFN'))
            self.output2.flush()

    def step(self, mcs):
        if OutputData:
            P = len(self.cell_list_by_type(self.I2))/self.shared_steppable_vars['InitialNumberCells']
            Time = mcs * hours_to_mcs
            ODEV = self.sbml.IFNModel['V']
            ODEH =  self.sbml.IFNModel['P']
            ODEP = self.sbml.IFNModel['P']
            ODEIFNe = self.sbml.IFNModel['IFNe'] * P
            ODESTATP = self.sbml.IFNModel['STATP']
            ODEIRF7 =  self.sbml.IFNModel['IRF7']
            ODEIRF7P =  self.sbml.IFNModel['IRF7P']
            ODEIFN =  self.sbml.IFNModel['IFN']
            self.output1.write("%f,%f,%f,%f,%f,%f,%f,%f,%f\n" %
                               (Time,ODEV,ODEH,ODEP,ODEIFNe,ODESTATP,ODEIRF7,ODEIRF7P,ODEIFN))
            self.output1.flush()

            L = len(self.cell_list_by_type(self.I2))
            CC3DP =  L / self.shared_steppable_vars['InitialNumberCells']
            CC3DV = 0.0
            CC3DH = 0.0
            CC3DSTATP = 0.0
            CC3DIRF7 = 0.0
            CC3DIRF7P = 0.0
            CC3DIFN = 0.0
            for cell in self.cell_list_by_type(self.I2):
                CC3DV += cell.sbml.VModel['V'] / L
                CC3DH += cell.sbml.VModel['H'] / L
                CC3DSTATP += cell.sbml.IModel['STATP'] / L
                CC3DIRF7 += cell.sbml.IModel['IRF7'] / L
                CC3DIRF7P += cell.sbml.IModel['IRF7P'] / L
                CC3DIFN += cell.sbml.IModel['IFN'] / L
            CC3DIFNe_Scalar = self.shared_steppable_vars['ExtracellularIFN_Scalar']\
                              / self.shared_steppable_vars['InitialNumberCells']
            CC3DIFNe_Field = self.shared_steppable_vars['ExtracellularIFN_Field'] \
                              / self.shared_steppable_vars['InitialNumberCells']
            self.output2.write("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" %
                               (Time,CC3DV,CC3DH,CC3DP,CC3DIFNe_Scalar,CC3DIFNe_Field,CC3DSTATP,CC3DIRF7,
                               CC3DIRF7P,CC3DIFN))
            self.output2.flush()