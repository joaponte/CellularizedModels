from cc3d.core.PySteppables import *
import numpy as np

plot_ODEModel = True
plot_CellModel = True
how_to_determine_IFNe = 1 # Determines the IFNe from the ODE model (1) from Cell model as scalar (2) or from field (3)

min_to_mcs = 1  # min/mcs
hours_to_mcs = min_to_mcs / 60.0 # hours/mcs
hours_to_simulate = 50.0

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

##Parameters
k11 = 0.0
k12 = 9.746
k13 = 12.511
k14 = 13.562
k21 = 10.385
t2 = 3.481
k61 = 0.635
k71 = 1.537
k72 = 47.883
k73 = 0.197
k31 = 45.922
k32 = 5.464
k33 = 0.068
t3 = 0.3
k41 = 0.115
k42 = 1.053
t4 = 0.75
k51 = 0.202
t5 = 0.3
n = 3.0
RIGI = 1.0

## Virus Model
def VModel(IFNe,H,V):
    dV = H*k71*V/(1.0+k72*IFNe*7E-5) - k73*V
    dH = - H*k61*V
    return(dV,dH)

## IFN Model
def IModel(IFNe,H,V, STATP,IRF7,IRF7P,IFN):
    dSTATP = H*k31*IFNe/(k32+k33*IFNe) - t3*STATP
    dIRF7 = H*(k41*STATP+k42*IRF7P) - t4*IRF7
    dIRF7P = H*k51*IRF7 - t5*IRF7P
    dIFN = H*(k11*RIGI*V+k12*(V**n)/(k13+(V**n))+k14*IRF7P) - k21*IFN
    return(dSTATP,dIRF7,dIRF7P,dIFN)

class ODEModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Load Original ODE Model
        self.add_free_floating_antimony(model_string=IFNModel_string, model_name='IFNModel',
                                        step_size=hours_to_mcs)

    def step(self,mcs):
        self.timestep_sbml()

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Store Initial Number of Cells
        self.InitialNumberCells = len(self.cell_list)

        # Store IFNe Scalar
        self.ExtracellularIFN_Scalar = 0.0

        # Set IFNe diffusion parameters
        self.get_xml_element('IFNe_dc').cdata = IFNe_diffusion_coefficient * min_to_mcs
        self.get_xml_element('IFNe_decay').cdata = t2 * hours_to_mcs

        # Set Max Simulation Steps
        self.get_xml_element('simulation_steps').cdata = hours_to_simulate / hours_to_mcs

        # Set secretors
        self.secretorIFN = self.get_field_secretor("IFNe")

        # Set Cell Properties
        for cell in self.cell_list:
            cell.dict['V'] = 6.9e-8
            cell.dict['H'] = 1.0
            cell.dict['STATP'] = 0.0
            cell.dict['IRF7'] = 0.0
            cell.dict['IRF7P'] = 0.0
            cell.dict['IFN'] = 0.0

    def step(self, mcs):
        ## Production of IFNe
        # E2b: IFN -> IFNe; k21 * IFN ;
        self.total_IFNe_production = 0.0
        k21C = k21 * hours_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            p = k21C * cell.dict['IFN']
            self.secretorIFN.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.total_IFNe_production += k21C * cell.dict['IFN']

        ## Decay of IFNe
        # E3a: IFNe -> ; t2*IFNe ;
        t2C = t2 * hours_to_mcs
        self.total_decay = t2C * self.ExtracellularIFN_Scalar
        ## Update Scalar IFNe
        self.ExtracellularIFN_Scalar += self.total_IFNe_production - self.total_decay

        ## P to D transition
        # E7a: P -> ; P * k61 * V;
        for cell in self.cell_list_by_type(self.I2):
            k61C = k61 * hours_to_mcs
            V = cell.dict['V']
            H = cell.dict['H']
            r = k61C * V * (1-H)
            p_I2toD = 1.0 - np.exp(-r)
            if np.random.random() < p_I2toD:
                cell.type = self.DEAD

        # Integrate Cell Quantities
        for cell in self.cell_list_by_type(self.I2):
            if how_to_determine_IFNe == 1:
                IFNe = self.sbml.IFNModel['IFNe']
            if how_to_determine_IFNe == 2:
                IFNe = self.ExtracellularIFN_Scalar / self.InitialNumberCells
            if how_to_determine_IFNe == 3:
                uptake_probability = 0.0000001
                uptake = self.secretorIFN.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
                IFNe = abs(uptake.tot_amount) / uptake_probability
                self.secretorIFN.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)
            V = cell.dict['V']
            H = cell.dict['H']
            STATP = cell.dict['STATP']
            IRF7 = cell.dict['IRF7']
            IRF7P = cell.dict['IRF7P']
            IFN = cell.dict['IFN']
            [dV,dH] = VModel(IFNe,H,V)
            [dSTATP, dIRF7,dIRF7P, dIFN] = IModel(IFNe,H,V, STATP,IRF7,IRF7P,IFN)
            cell.dict['V'] += dV * hours_to_mcs
            cell.dict['H'] += dH * hours_to_mcs
            cell.dict['STATP'] += dSTATP * hours_to_mcs
            cell.dict['IRF7'] += dIRF7 * hours_to_mcs
            cell.dict['IRF7P'] += dIRF7P * hours_to_mcs
            cell.dict['IFN'] += dIFN * hours_to_mcs

class IFNPlotSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Store Initial Number of Cells
        self.InitialNumberCells = len(self.cell_list)

        # Measure IFNe Scalar
        self.ExtracellularIFN_Scalar = 0.0

        # Set secretors
        self.secretorIFN = self.get_field_secretor("IFNe")

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
            P =  len(self.cell_list_by_type(self.I2))/self.InitialNumberCells
            self.plot_win1.add_data_point("ODEV", mcs * hours_to_mcs, self.sbml.IFNModel['V'])
            self.plot_win2.add_data_point("ODEH", mcs * hours_to_mcs, self.sbml.IFNModel['P'])
            self.plot_win3.add_data_point("ODEP", mcs * hours_to_mcs, self.sbml.IFNModel['P'])
            self.plot_win4.add_data_point("ODEIFNe", mcs * hours_to_mcs, self.sbml.IFNModel['IFNe'] * P)
            self.plot_win5.add_data_point("ODESTATP", mcs * hours_to_mcs, self.sbml.IFNModel['STATP'])
            self.plot_win6.add_data_point("ODEIRF7", mcs * hours_to_mcs, self.sbml.IFNModel['IRF7'])
            self.plot_win7.add_data_point("ODEIRF7P", mcs * hours_to_mcs, self.sbml.IFNModel['IRF7P'])
            self.plot_win8.add_data_point("ODEIFN", mcs * hours_to_mcs, self.sbml.IFNModel['IFN'])

        if plot_CellModel:
            ## Production of IFNe
            # E2b: IFN -> IFNe; k21 * IFN ;
            self.total_IFNe_production = 0.0
            k21C = k21 * hours_to_mcs
            for cell in self.cell_list_by_type(self.I2):
                intracellularIFN = cell.dict['IFN']
                self.total_IFNe_production += k21C * intracellularIFN

            ## Decay of IFNe
            # E3a: IFNe -> ; t2*IFNe ;
            I = self.ExtracellularIFN_Scalar
            t2C = t2 * hours_to_mcs
            self.total_decay = t2C * I
            ## Update Scalar IFNe
            self.ExtracellularIFN_Scalar += self.total_IFNe_production - self.total_decay

            self.ExtracellularIFN_Field = 0.0
            for cell in self.cell_list:
                uptake_probability = 0.0000001
                uptake = self.secretorIFN.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
                self.ExtracellularIFN_Field += abs(uptake.tot_amount) / uptake_probability
                self.secretorIFN.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

            L = len(self.cell_list_by_type(self.I2))
            avgV = 0.0
            avgH = 0.0
            avgSTATP = 0.0
            avgIRF7 = 0.0
            avgIRF7P = 0.0
            avgIFN = 0.0
            for cell in self.cell_list_by_type(self.I2):
                avgV += cell.dict['V'] / L
                avgH += cell.dict['H'] / L
                avgSTATP += cell.dict['STATP'] / L
                avgIRF7 += cell.dict['IRF7'] / L
                avgIRF7P += cell.dict['IRF7P'] / L
                avgIFN += cell.dict['IFN']/ L

            self.plot_win1.add_data_point("CC3DV", mcs * hours_to_mcs, avgV)
            self.plot_win2.add_data_point("CC3DH", mcs * hours_to_mcs, avgH)
            self.plot_win3.add_data_point("CC3DP", mcs * hours_to_mcs,
                                          L / self.InitialNumberCells)
            self.plot_win4.add_data_point("CC3DIFNe_Scalar", mcs * hours_to_mcs,
                                          self.ExtracellularIFN_Scalar /self.InitialNumberCells)
            self.plot_win4.add_data_point("CC3DIFNe_Field", mcs * hours_to_mcs,
                                          self.ExtracellularIFN_Field/self.InitialNumberCells)
            self.plot_win5.add_data_point("CC3DSTATP", mcs * hours_to_mcs, avgSTATP)
            self.plot_win6.add_data_point("CC3DIRF7", mcs * hours_to_mcs, avgIRF7)
            self.plot_win7.add_data_point("CC3DIRF7P", mcs * hours_to_mcs, avgIRF7P)
            self.plot_win8.add_data_point("CC3DIFN", mcs * hours_to_mcs, avgIFN)