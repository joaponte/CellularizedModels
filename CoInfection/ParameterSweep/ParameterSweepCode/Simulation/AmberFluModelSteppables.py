from cc3d.core.PySteppables import *
import numpy as np
import os
import Parameters

Data_writeout_ODEs = False
Data_writeout_CellularModel = True

min_to_mcs = 10.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 10.0 #10 in the original model

# production_multiplier = Parameters.P
# diffusion_multiplier = Parameters.D
replicate = Parameters.R
infection_time = Parameters.T
K_V = Parameters.K

beta_ode = 2.4 * 10**(-4)
p_ode = 1.6
c_ode = 13.0
k_ode = 4.0
delta_d_ode = 1.6 * 10**6
K_delta_ode = 4.5 * 10**5
T0_ode = 1.0 * 10**7

class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Uptading max simulation steps using scaling factor to simulate 10 days
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs

        # set initial model parameters
        self.initial_uninfected = len(self.cell_list_by_type(self.U))

        self.get_xml_element('virus_dc').cdata = 1.0 #* diffusion_multiplier
        self.get_xml_element('virus_decay').cdata = c_ode * days_to_mcs

        self.get_xml_element('virusB_dc').cdata = 1.0
        self.get_xml_element('virusB_decay').cdata = c_ode * days_to_mcs

    def step(self, mcs):
        start_time = infection_time / days_to_mcs
        if mcs == int(start_time):
            virus_field_amount = 75.0 / (self.dim.x * self.dim.y)
            for x in range(0,self.dim.x):
                for y in range(0,self.dim.y):
                    self.field.VirusB[x, y, 0] = virus_field_amount

        # Measure total extracellular VirusA
        secretorA = self.get_field_secretor("Virus")
        self.Virus_Field = 0
        for cell in self.cell_list:
            uptake_probability = 0.0000001
            uptake = secretorA.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            V = abs(uptake.tot_amount) / uptake_probability
            self.Virus_Field += V
            secretorA.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)
            
        # Measure total extracellular VirusB
        secretorB = self.get_field_secretor("VirusB")
        self.Virus_FieldB = 0
        for cell in self.cell_list:
            uptake_probability = 0.0000001
            uptake = secretorB.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            V = abs(uptake.tot_amount) / uptake_probability
            self.Virus_FieldB += V
            secretorB.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

        # Transition rule from U to I1
        for cell in self.cell_list_by_type(self.U):
            # Determine V from the virus field
            b = beta_ode * self.initial_uninfected * days_to_mcs
            uptake_probability = 0.0000001

            uptake = secretorA.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            VA = abs(uptake.tot_amount) / uptake_probability
            secretorA.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

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
        k = k_ode * days_to_mcs
        p_T1oI2 = k
        for cell in self.cell_list_by_type(self.I1):
            if np.random.random() < p_T1oI2:
                cell.type = self.I2

        # Transition rule from I1B to I2B
        k = k_ode * days_to_mcs
        p_T1BoI2B = k
        for cell in self.cell_list_by_type(self.I1B):
            if np.random.random() < p_T1BoI2B:
                cell.type = self.I2B

        I2 = len(self.cell_list_by_type(self.I2))
        I2B = len(self.cell_list_by_type(self.I2B))
        # Transition rule from I2 to D
        K_delta = K_delta_ode / T0_ode * self.initial_uninfected
        delta_d = delta_d_ode / T0_ode * self.initial_uninfected
        # VirusB infects CD8 for VirusA
        # effect = K_V / (K_V + self.Virus_FieldB)
        effect = 1.0
        p_T2toD = effect * delta_d / (K_delta + I2) * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            if np.random.random() < p_T2toD:
                cell.type = self.DEAD

        # Transition rule from I2B to DB
        K_delta = K_delta_ode / T0_ode * self.initial_uninfected
        delta_d = delta_d_ode / T0_ode * self.initial_uninfected
        # VirusA infects CD8 for VirusB
        effect = K_V / (K_V + self.Virus_Field)
        p_T2BtoDB = effect * delta_d / (K_delta + I2B) * days_to_mcs
        for cell in self.cell_list_by_type(self.I2B):
            if np.random.random() < p_T2BtoDB:
                cell.type = self.DEADB

        # Production of extracellular virus A
        p = p_ode / self.initial_uninfected * T0_ode * days_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            secretorA.secreteInsideCellTotalCount(cell, p / cell.volume)

        # Production of extracellular virus A
        p = p_ode / self.initial_uninfected * T0_ode * days_to_mcs
        for cell in self.cell_list_by_type(self.I2B):
            secretorB.secreteInsideCellTotalCount(cell, p / cell.volume)

class Data_OutputSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        if Data_writeout_CellularModel:
            folder_path = '/Users/Josua/Downloads/AmberFluModelv3/'
            # folder_path = '/N/u/joaponte/Carbonate/FluModel/Output/'
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)

            file_name2 = 'cellularizedmodel_%.2f_%i_%i.txt' % (infection_time, K_V, replicate)
            self.output2 = open(folder_path + file_name2, 'w')
            self.output2.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                'Time', 'U', 'AI1', 'AI2', 'AD', 'BI1', 'BI2', 'BD', 'VA', 'VB'))
            self.output2.flush()

    def step(self, mcs):

        if Data_writeout_CellularModel:
            # Record variables from Cellularized Model
            d = mcs * days_to_mcs
            U = len(self.cell_list_by_type(self.U))
            I1 = len(self.cell_list_by_type(self.I1))
            I2 = len(self.cell_list_by_type(self.I2))
            DA = len(self.cell_list_by_type(self.DEAD))
            I1B = len(self.cell_list_by_type(self.I1B))
            I2B = len(self.cell_list_by_type(self.I2B))
            DB = len(self.cell_list_by_type(self.DEADB))

            secretorA = self.get_field_secretor("Virus")
            secretorB = self.get_field_secretor("VirusB")

            self.Virus_Field = 0
            for cell in self.cell_list:
                uptake_probability = 0.0000001
                uptake = secretorA.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
                V = abs(uptake.tot_amount) / uptake_probability
                self.Virus_Field += V
                secretorA.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

            self.Virus_FieldB = 0
            for cell in self.cell_list:
                uptake_probability = 0.0000001
                uptake = secretorB.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
                V = abs(uptake.tot_amount) / uptake_probability
                self.Virus_FieldB += V
                secretorB.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

            self.output2.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (
            d, U, I1, I2, DA, I1B, I2B, DB, self.Virus_Field, self.Virus_FieldB))
            self.output2.flush()

    def finish(self):
        if Data_writeout_CellularModel:
            self.output2.close()
