from cc3d.core.PySteppables import *
import numpy as np
import os
import Parameters

min_to_mcs = 10.0  # min/mcs
hours_to_mcs = min_to_mcs / 60.0 # hours/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
hours_to_simulate = 50.0  # 10 in the original model

Replicate = Parameters.R
Multiplier = Parameters.M

def Virus_Model(H,V,IFN):
    k61 = 0.635
    k71 = 1.537
    k72 = 47.883
    k73 = 0.197
    dV = H*(k71*V)/(1.0+k72*IFN*7E-5) - k73*V
    dH = -H*k61*V
    return(dH,dV)

def IFN_Model(H,V, IFNe,STATP,IRF7,IRF7P,IFN):
    k12 = 9.746
    k13 = 12.511
    k14 = 13.562
    k21 = 10.385
    k31 = 45.922
    k32 = 5.464
    k33 = 0.068
    t3  = 0.3
    k41 = 0.115
    k42 = 1.053
    t4  = 0.75
    k51 = 0.202
    t5  = 0.3
    n = 3.0
    dSTATP = H*k31*IFNe/(k32+k33*IFNe)-t3*STATP
    dIRF7 = H*(k41*STATP+k42*IRF7P)-t4*IRF7
    dIRF7P = H*k51*IRF7-t5*IRF7P
    dIFN = H*(k12*(V**n)/(k13+(V**n))+k14*IRF7P) - k21*IFN
    return(dSTATP,dIRF7,dIRF7P,dIFN)

class ModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Flu Model Parameters
        self.beta = 2.4 * 10 ** (-4)
        self.k = 4.0
        self.c = 13.0

        # IFN Model Parameters
        self.k61 = 0.635
        self.k21 = 10.385
        self.t2 = 3.481
        self.k73 = 0.197

        # Diffusion Parameters
        self.get_xml_element('simulation_steps').cdata = hours_to_simulate / hours_to_mcs
        self.get_xml_element('IFNe_decay').cdata = self.t2 * hours_to_mcs
        self.get_xml_element('virus_decay').cdata = self.c * days_to_mcs

        # Tracking Variables
        self.InitialNumberCells = len(self.cell_list)

        for cell in self.cell_list:
            cell.dict['H'] = 1.0
            cell.dict['V'] = 0.0
            cell.dict['STATP'] = 0.0
            cell.dict['IRF7'] = 0.0
            cell.dict['IRF7P'] = 0.0
            cell.dict['IFN'] = 0.0

        cell = self.cell_field[self.dim.x // 2, self.dim.y // 2, 0]
        cell.type = self.I1
        cell.dict['V'] = 6.9e-8

    def step(self, mcs):
        secretorIFN = self.get_field_secretor("IFNe")
        secretorV = self.get_field_secretor("Virus")

        # Update Internal Cellular Models
        for cell in self.cell_list:
            STATP = cell.dict['STATP']
            IRF7 = cell.dict['IRF7']
            IRF7P = cell.dict['IRF7P']
            IFN = cell.dict['IRF7P']
            H = cell.dict['H']
            V = cell.dict['V']
            uptake_probability = 0.0000001
            uptake = secretorIFN.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            IFNe = abs(uptake.tot_amount) / uptake_probability
            secretorIFN.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)
            (dSTATP,dIRF7,dIRF7P,dIFN) = IFN_Model(H,V,IFNe,STATP,IRF7,IRF7P,IFN)
            (dH, dV) = Virus_Model(H,V,IFN)
            cell.dict['STATP'] = dSTATP * hours_to_mcs + STATP
            cell.dict['IRF7'] = dIRF7 * hours_to_mcs + IRF7
            cell.dict['IRF7P'] = dIRF7P * hours_to_mcs + IRF7P
            cell.dict['IFN'] = dIFN * hours_to_mcs + IFN
            cell.dict['H'] = dH * hours_to_mcs + H
            cell.dict['V'] = dV * hours_to_mcs + V

        # Transition from T to I1
        for cell in self.cell_list_by_type(self.U):
            # Determine V from the virus field
            b = self.beta * self.InitialNumberCells * days_to_mcs
            uptake_probability = 0.0000001
            uptake = secretorV.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            V = abs(uptake.tot_amount) / uptake_probability
            secretorV.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)
            p_UtoI1 = b * V
            if np.random.random() < p_UtoI1:
                cell.type = self.I1
                cell.dict['V'] = 6.9e-8

        # Transition from I1 to I2
        for cell in self.cell_list_by_type(self.I1):
            k = self.k * days_to_mcs
            p_T1oI2 = k
            if np.random.random() < p_T1oI2:
                cell.type = self.I2

        # Transition from I2 to Dead
        for cell in self.cell_list_by_type(self.I2):
            k61 = self.k61 * hours_to_mcs
            V = cell.dict['V']
            H = cell.dict['H']
            p_I2toD = k61 * V * (1-H)
            if np.random.random() < p_I2toD:
                cell.type = self.DEAD

        ## Production of extracellular IFN
        k21 = self.k21 * hours_to_mcs
        for cell in self.cell_list_by_type(self.U,self.I1,self.I2):
            secretorIFN.secreteInsideCellTotalCount(cell, k21 * cell.dict['IFN'] / cell.volume)

        ## Production of extracellular virus
        k73 = self.k73 * hours_to_mcs
        for cell in self.cell_list_by_type(self.I2):
            secretorV.secreteInsideCellTotalCount(cell, k73 * cell.dict['V'] * 1094460.28 / cell.volume)

class PlaqueAssaySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        folder_path = '/Users/Josua/Data/'
        # folder_path = '/N/u/joaponte/Carbonate/FluModel/Output/'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_name = 'PlaqueAssay_%s_%.3f_%i.txt' % ('k61',Multiplier,Replicate)
        self.output = open(folder_path + file_name, 'w')
        self.output.write("%s,%s,%s,%s,%s\n" % ('Time','avgI1rd', 'avgI2rd', 'avgDrd', 'Beff'))
        self.output.flush()

        # Parameters for Measuring Effective Infectivity
        self.previousT = 0.0

    def step(self, mcs):
        h = mcs * hours_to_mcs

        avgI1rd = 0.0
        num_I1 = len(self.cell_list_by_type(self.I1))
        for cell in self.cell_list_by_type(self.I1):
            xCOM = cell.xCOM
            yCOM = cell.yCOM
            avgI1rd += sqrt((self.dim.x/2.0 - xCOM)**2 + (self.dim.y/2.0-yCOM)**2) / num_I1

        avgI2rd = 0.0
        num_I2 = len(self.cell_list_by_type(self.I2))
        for cell in self.cell_list_by_type(self.I2):
            xCOM = cell.xCOM
            yCOM = cell.yCOM
            avgI2rd += sqrt((self.dim.x/2.0 - xCOM)**2 + (self.dim.y/2.0-yCOM)**2) / num_I2

        avgDrd = 0.0
        num_D = len(self.cell_list_by_type(self.DEAD))
        for cell in self.cell_list_by_type(self.DEAD):
            xCOM = cell.xCOM
            yCOM = cell.yCOM
            avgDrd += sqrt((self.dim.x/2.0 - xCOM)**2 + (self.dim.y/2.0-yCOM)**2) / num_D

        secretorV = self.get_field_secretor("Virus")
        Beff = 0.0
        num_T = len(self.cell_list_by_type(self.U))
        dT = abs(num_T - self.previousT)
        self.previousT = num_T
        self.ExtracellularVirus_Field = 0.0
        for cell in self.cell_list:
            uptake_probability = 0.0000001
            uptake = secretorV.uptakeInsideCellTotalCount(cell, 1E6, uptake_probability)
            self.ExtracellularVirus_Field  += abs(uptake.tot_amount) / uptake_probability
            secretorV.secreteInsideCellTotalCount(cell, abs(uptake.tot_amount) / cell.volume)

        if self.ExtracellularVirus_Field:
            Beff = dT / (num_T*self.ExtracellularVirus_Field*hours_to_mcs)

        self.output.write("%.2f,%.2f,%.2f,%.2f,%f\n" % (h, avgI1rd , avgI2rd , avgDrd, Beff))
        self.output.flush()
