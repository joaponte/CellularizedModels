from cc3d.core.PySteppables import *
import numpy as np
import os

plot_SingleVirusODEModel = False
plot_CoinfectionODEModel = False
plot_CoinfectionCellModel = False
Data_writeout = False

## How to determine V
# 1 pulls from the scalar virus from the ODE original model (no feedback in the cellular model)
# 2 pulls from the scalar virus from the cellular model (feedback in the cellular model but no field)
# 3 pulls from the virus field
how_to_determine_V = 3

min_to_mcs = 10.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 4.0  # 10 in the original model

'''Smith AP, Moquin DJ, Bernhauerova V, Smith AM. Influenza virus infection model with density dependence 
supports biphasic viral decay. Frontiers in microbiology. 2018 Jul 10;9:1554.'''

FluModelString = '''        
        model ambersmithsimple()
        
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

CoinfectionModelString = '''
        model coinfection()
         //State Variables and Transitions for Virus A
        V1: -> T   ; -beta * VA * T ;                               // Susceptible Cells
        V2: -> I1A ;  beta * VA * T - k * I1A ;                     // Early Infected Cells
        V3: -> I2A ;  k * I1A - delta_d * I2A / (K_delta + I2A) ;   // Late Infected Cells
        V4: -> VA  ;  p * I2A - c * VA;                             // Extracellular Virus A
        V5: -> DA  ;  delta_d * I2A / (K_delta + I2A) ;             // Dead Cells

         //State Variables and Transitions for Virus B
        V6: -> T   ; -beta * VB * T ;                               // Susceptible Cells
        V7: -> I1B ;  beta * VB * T - k * I1B ;                     // Early Infected Cells
        V8: -> I2B ;  k * I1B - delta_d * I2B / (K_delta + I2B) ;   // Late Infected Cells
        V9: -> VB  ;  p * I2B - c * VB;                             // Extracellular Virus B
        V10: -> DB  ;  delta_d * I2B / (K_delta + I2B) ;             // Dead Cells

        //Parameters
        beta = 2.4* 10^(-4) ;                                       // Virus Infective
        p = 1.6 ;                                                   // Virus Production
        c = 13.0 ;                                                  // Virus Clearance
        k = 4.0 ;                                                   // Eclipse phase
        delta_d = 1.6 * 10^6 ;                                      // Infected Cell Clearance
        K_delta = 4.5 * 10^5 ;                                      // Half Saturation Constant   

        // Initial Conditions ;
        T0 = 1.0*10^7;
        T = T0  ;                                                   // Initial Number of Uninfected Cells
        I1 = 75.0 ;                                                 // Initial Number of Infected Cells
end'''

class AmberFluModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Uptading max simulation steps using scaling factor to simulate 10 days
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs

        # Adding free floating single antimony model
        self.add_free_floating_antimony(model_string=FluModelString, model_name='FluModel',
                                        step_size=days_to_mcs)

        # Changing initial values according to discussions with Amber Smith
        self.sbml.FluModel['I1'] = 0.0
        self.sbml.FluModel['V'] = 75.0

        # Adding free floating coinfection antimony model
        self.add_free_floating_antimony(model_string=CoinfectionModelString, model_name='CoinfectionModel',
                                        step_size=days_to_mcs)

        # Changing initial values according to discussions with Amber Smith
        self.sbml.CoinfectionModel['I1'] = 0.0
        self.sbml.CoinfectionModel['VA'] = 75.0
        self.sbml.CoinfectionModel['VB'] = 75.0

    def step(self, mcs):
        self.timestep_sbml()

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # set initial model parameters
        self.initial_uninfected = len(self.cell_list_by_type(self.U))
        self.ExtracellularVirusA = self.sbml.CoinfectionModel['VA']
        self.get_xml_element('virusA_decay').cdata = self.sbml.CoinfectionModel['c'] * days_to_mcs
        self.ExtracellularVirusB = self.sbml.coinfection['VB']
        self.get_xml_element('virusB_decay').cdata = self.sbml.CoinfectionModel['c'] * days_to_mcs

    def step(self, mcs):
        # Transition rule from U to I1
        secretorA = self.get_field_secretor("VirusA")
        secretorB = self.get_field_secretor("VirusB")
        for cell in self.cell_list_by_type(self.U):
            # Determine V from scalar virus from the ODE
            if how_to_determine_V == -1:
                b = self.sbml.CoinfectionModel['beta'] * self.sbml.CoinfectionModel['T0'] * days_to_mcs
                VA = self.sbml.CoinfectionModel['VA'] / self.sbml.CoinfectionModel['T0']
                VB = self.sbml.CoinfectionModel['VB'] / self.sbml.CoinfectionModel['T0']

            # Determine V from scalar virus from the cellular model
            if how_to_determine_V == 0:
                b = self.sbml.CoinfectionModel['beta'] * self.initial_uninfected * days_to_mcs
                VA = self.ExtracellularVirusA / self.initial_uninfected
                VB = self.ExtracellularVirusB / self.initial_uninfected

            # Determine V from the virus field
            if how_to_determine_V == 1:
                b = self.sbml.CoinfectionModel['beta'] * self.initial_uninfected * days_to_mcs
                uptake_probability = 0.0000001

                # Measure Virus A
                uptake = secretorA.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
                VA = abs(uptake.tot_amount) / uptake_probability
                secretorA.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

                # Measure Virus B
                uptakeB = secretorB.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
                VB = abs(uptakeB.tot_amount) / uptake_probability
                secretorB.secreteInsideCellTotalCount(cell, abs(uptakeB.tot_amount) / cell.volume)

            # Calculate the probability of infection of individual cells based on the amount of virus PER cell
            # Transition rule from T to I2
            p_UtoI1A = b * VA
            if np.random.random() < p_UtoI1A:
                cell.type = self.I1

            # Transition rule from T to I2B
            p_UtoI1A = b * VB
            if np.random.random() < p_UtoI1A:
                cell.type = self.I1B

        # Transition rule from I1 to I2
        k = self.sbml.CoinfectionModel['k'] * days_to_mcs
        p_T1oI2 = k
        for cell in self.cell_list_by_type(self.I1):
            if np.random.random() < p_T1oI2:
                cell.type = self.I2

        # Transition rule from I1B to I2B
        k = self.sbml.CoinfectionModel['k'] * days_to_mcs
        p_T1BoI2B = k
        for cell in self.cell_list_by_type(self.I1B):
            if np.random.random() < p_T1BoI2B:
                cell.type = self.I2B

        # Transition rule from I2 to D
        K_delta = self.sbml.CoinfectionModel['K_delta'] / self.sbml.CoinfectionModel['T0'] * self.initial_uninfected
        delta_d = self.sbml.CoinfectionModel['delta_d'] / self.sbml.CoinfectionModel['T0'] * self.initial_uninfected
        I2 = len(self.cell_list_by_type(self.I2))
        p_T2toD = delta_d / (K_delta + I2) * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            if np.random.random() < p_T2toD:
                cell.type = self.DEAD

        # Transition rule from I2B to DB
        K_delta = self.sbml.CoinfectionModel['K_delta'] / self.sbml.CoinfectionModel['T0'] * self.initial_uninfected
        delta_d = self.sbml.CoinfectionModel['delta_d'] / self.sbml.CoinfectionModel['T0'] * self.initial_uninfected
        I2B = len(self.cell_list_by_type(self.I2B))
        p_T2BtoDB = delta_d / (K_delta + I2B) * days_to_mcs
        for cell in self.cell_list_by_type(self.I2B):
            if np.random.random() < p_T2BtoDB:
                cell.type = self.DEADB

        # Production of extracellular virus A
        secretor = self.get_field_secretor("VirusA")
        VA = self.ExtracellularVirusA
        p = self.sbml.CoinfectionModel['p'] / self.initial_uninfected * self.sbml.CoinfectionModel['T0'] * days_to_mcs
        c = self.sbml.CoinfectionModel['c'] * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            release = secretor.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.ExtracellularVirusA += release.tot_amount
        self.ExtracellularVirusA -= c * VA
        self.shared_steppable_vars['ScalarVirusA'] = self.ExtracellularVirusA

        # Production of extracellular virus B
        secretor = self.get_field_secretor("VirusB")
        VB = self.ExtracellularVirusB
        p = self.sbml.CoinfectionModel['p'] / self.initial_uninfected * self.sbml.CoinfectionModel['T0'] * days_to_mcs
        c = self.sbml.CoinfectionModel['c'] * days_to_mcs
        for cell in self.cell_list_by_type(self.I2B):
            release = secretor.secreteInsideCellTotalCount(cell, p / cell.volume)
            self.ExtracellularVirusB += release.tot_amount
        self.ExtracellularVirusB -= c * VB
        self.shared_steppable_vars['ScalarVirusB'] = self.ExtracellularVirusB

        # Measure extracellular virus field A
        self.Plot_ExtracellularVirus_Field = 0
        for cell in self.cell_list:
            uptake_probability = 0.0000001
            uptake = secretor.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            V = abs(uptake.tot_amount) / uptake_probability
            self.Plot_ExtracellularVirus_Field += V
            secretorA.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)
        self.shared_steppable_vars['FieldVirusA'] = self.Plot_ExtracellularVirus_FieldA

        # Measure amount of extracellular virus field B
        self.Plot_ExtracellularVirus_FieldB = 0
        for cell in self.cell_list:
            uptake_probability = 0.0000001
            uptake = secretor.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            V = abs(uptake.tot_amount) / uptake_probability
            self.Plot_ExtracellularVirus_FieldB += V
            secretorB.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)
        self.shared_steppable_vars['FieldVirusB'] = self.Plot_ExtracellularVirus_FieldB

class Data_OutputSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        if Data_writeout:
            folder_path = '/Users/Josua/AmberFluModelData/'
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)

            file_name = 'CoinfectionODE_output.txt'
            self.output = open(folder_path + file_name, 'w')
            self.output.write(
                "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ('Time', 'T', 'I1A', 'I2A', 'DA', 'VA', 'I1B', 'I2B', 'DB', 'VB'))
            self.output.flush()

            file_name = 'CoinfectionCellular_output.txt'
            self.output2 = open(folder_path + file_name, 'w')
            self.output2.write(
                "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ('Time', 'T', 'I1A', 'I2A', 'DA', 'VA', 'I1B', 'I2B', 'DB', 'VB'))
            self.output2.flush()

    def step(self, mcs):
        if Data_writeout:
            # Record variables from ODE model
            d = mcs * days_to_mcs
            AT = self.sbml.coinfection['T']
            AI1A = self.sbml.coinfection['I1A']
            AI2A = self.sbml.coinfection['I2A']
            ADA = self.sbml.coinfection['DA']
            AVA = self.sbml.coinfection['VA']
            AI1B = self.sbml.coinfection['I1B']
            AI2B = self.sbml.coinfection['I2B']
            ADB = self.sbml.coinfection['DB']
            AVB = self.sbml.coinfection['VB']

            self.output.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (
                d, AT, AI1A, AI2A, ADA, AVA, AI1B, AI2B, ADB, AVB))
            self.output.flush()

            # Record variables from Cellularized Model
            d = mcs * days_to_mcs
            CU = len(self.cell_list_by_type(self.U))
            CI1A = len(self.cell_list_by_type(self.I1))
            CI2A = len(self.cell_list_by_type(self.I2))
            CDA = len(self.cell_list_by_type(self.DEAD))
            CI1B = len(self.cell_list_by_type(self.I1B))
            CI2B = len(self.cell_list_by_type(self.I2B))
            CDB = len(self.cell_list_by_type(self.DEADB))

            if how_to_determine_V == 3:
                CVA = self.shared_steppable_vars['FieldVirusA']
                CVB = self.shared_steppable_vars['FieldVirusB']
            else:
                CVA = self.shared_steppable_vars['ScalarVirusA']
                CVB = self.shared_steppable_vars['ScalarVirusB']

            self.output.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (
                d, CU, CI1A, CI2A, CDA, CVA, CI1B, CI2B, CDB, CVB))
            self.output.flush()

    def finish(self):
        if Data_writeout:
            self.output.close()

class PlotSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Initialize Graphic Window for Amber Smith ODE model
        self.initial_uninfected = len(self.cell_list_by_type(self.U))

        if (plot_SingleVirusODEModel == True) or (plot_CoinfectionODEModel == True) or (plot_CoinfectionCellModel  == True):
            self.plot_win = self.add_new_plot_window(title='Cells',
                                                     x_axis_title='Days',
                                                     y_axis_title='Variables', x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=False)#, config_options={'legend': True})

            self.plot_win2 = self.add_new_plot_window(title='Virus',
                                                      x_axis_title='Days',
                                                      y_axis_title='Variables', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False)# ,config_options={'legend': True})

        if plot_SingleVirusODEModel == True:
            self.plot_win.add_plot("SiT", style='Dots', color='blue', size=5)
            self.plot_win.add_plot("SiI1", style='Dots', color='yellow', size=5)
            self.plot_win.add_plot("SiI2", style='Dots', color='red', size=5)
            self.plot_win.add_plot("SiD", style='Dots', color='purple', size=5)
            self.plot_win2.add_plot("SiV", style='Dots', color='blue', size=5)

        # Initialize Graphic Window for Coinfection ODE model
        if plot_CoinfectionODEModel == True:
            self.plot_win.add_plot("CoODET", style='Dots', color='blue', size=5)
            self.plot_win.add_plot("CoODEI1A", style='Dots', color='yellow', size=5)
            self.plot_win.add_plot("CoODEI2A", style='Dots', color='red', size=5)
            self.plot_win.add_plot("CoODEDA", style='Dots', color='purple', size=5)
            self.plot_win2.add_plot("CoODEVA", style='Dots', color='blue', size=5)
            self.plot_win.add_plot("CoODEI1B", style='Dots', color='green', size=5)
            self.plot_win.add_plot("CoODEI2B", style='Dots', color='magenta', size=5)
            self.plot_win.add_plot("CoODEDB", style='Dots', color='brown', size=5)
            self.plot_win2.add_plot("CoODEVB", style='Dots', color='cyan', size=5)

        if plot_CoinfectionCellModel == True:
            self.plot_win.add_plot("CoCC3DT", style='Lines', color='blue', size=5)
            self.plot_win.add_plot("CoCC3DI1A", style='Lines', color='yellow', size=5)
            self.plot_win.add_plot("CoCC3DI2A", style='Lines', color='red', size=5)
            self.plot_win.add_plot("CoCC3DDA", style='Lines', color='purple', size=5)
            self.plot_win2.add_plot("CoCC3DVA", style='Lines', color='blue', size=5)
            self.plot_win.add_plot("CoCC3DI1B", style='Lines', color='green', size=5)
            self.plot_win.add_plot("CoCC3DI2B", style='Lines', color='magenta', size=5)
            self.plot_win.add_plot("CoCC3DDB", style='Lines', color='brown', size=5)
            self.plot_win2.add_plot("CoCC3DVB", style='Lines', color='cyan', size=5)

    def step(self, mcs):
        if plot_SingleVirusODEModel == True:
            self.plot_win.add_data_point("SiT", mcs * days_to_mcs,
                                         self.sbml.FluModel['T'] / self.sbml.FluModel['T0'])
            self.plot_win.add_data_point("SiI1", mcs * days_to_mcs,
                                         self.sbml.FluModel['I1'] / self.sbml.FluModel['T0'])
            self.plot_win.add_data_point("SiI2", mcs * days_to_mcs,
                                         self.sbml.FluModel['I2'] / self.sbml.FluModel['T0'])
            self.plot_win.add_data_point("SiD", mcs * days_to_mcs,
                                         self.sbml.FluModel['D'] / self.sbml.FluModel['T0'])
            self.plot_win2.add_data_point("SiV", mcs * days_to_mcs, np.log10(self.sbml.FluModel['V']))

        if plot_CoinfectionODEModel == True:
            self.plot_win.add_data_point("CoODET", mcs * days_to_mcs,
                                          self.sbml.CoinfectionModel['T'] / self.sbml.CoinfectionModel['T0'])
            self.plot_win.add_data_point("CoODEI1A", mcs * days_to_mcs,
                                          self.sbml.CoinfectionModel['I1A'] / self.sbml.CoinfectionModel['T0'])
            self.plot_win.add_data_point("CoODEI2A", mcs * days_to_mcs,
                                          self.sbml.CoinfectionModel['I2A'] / self.sbml.CoinfectionModel['T0'])
            self.plot_win.add_data_point("CoODEDA", mcs * days_to_mcs,
                                          self.sbml.CoinfectionModel['DA'] / self.sbml.CoinfectionModel['T0'])
            self.plot_win2.add_data_point("CoODEVA", mcs * days_to_mcs,
                                          np.log10(self.sbml.CoinfectionModel['VA']))
            self.plot_win.add_data_point("CoODEI1B", mcs * days_to_mcs,
                                          self.sbml.coinfection['I1B'] / self.sbml.CoinfectionModel['T0'])
            self.plot_win.add_data_point("CoODEI2B", mcs * days_to_mcs,
                                          self.sbml.coinfection['I2B'] / self.sbml.CoinfectionModel['T0'])
            self.plot_win.add_data_point("CoODEDB", mcs * days_to_mcs,
                                          self.sbml.coinfection['DB'] / self.sbml.CoinfectionModel['T0'])
            self.plot_win2.add_data_point("CoODEVB", mcs * days_to_mcs, np.log10(self.sbml.CoinfectionModel['VB']))

        if plot_CoinfectionCellModel == True:
            self.Plot_ExtracellularVirusA = self.shared_steppable_vars['ScalarVirusA']
            self.Plot_ExtracellularVirusB = self.shared_steppable_vars['ScalarVirusB']
            self.Plot_ExtracellularVirus_FieldA = self.shared_steppable_vars['FieldVirusA']
            self.Plot_ExtracellularVirus_FieldB = self.shared_steppable_vars['FieldVirusB']

            self.plot_win.add_data_point("CoCC3DT", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.U)) / self.initial_uninfected)
            self.plot_win.add_data_point("CoCC3I1A", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.I1)) / self.initial_uninfected)
            self.plot_win.add_data_point("CoCC3I2A", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.I2)) / self.initial_uninfected)
            self.plot_win.add_data_point("CoCC3DA", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.DEAD)) / self.initial_uninfected)
            self.plot_win.add_data_point("CoCC3I1B", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.I1B)) / self.initial_uninfected)
            self.plot_win.add_data_point("CoCC3I2B", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.I2B)) / self.initial_uninfected)
            self.plot_win.add_data_point("CoCC3DB", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.DEADB)) / self.initial_uninfected)

            if how_to_determine_V == 3:
                self.plot_win4.add_data_point("CoCC3VA", mcs * days_to_mcs,
                                              np.log10(self.Plot_ExtracellularVirus_FieldA))
                self.plot_win4.add_data_point("CoCC3VB", mcs * days_to_mcs,
                                              np.log10(self.Plot_ExtracellularVirus_FieldB))
            else:
                self.plot_win4.add_data_point("CoCC3VA", mcs * days_to_mcs,
                                              np.log10(self.Plot_ExtracellularVirusA))
                self.plot_win4.add_data_point("CoCC3VB", mcs * days_to_mcs,
                                              np.log10(self.Plot_ExtracellularVirusB))