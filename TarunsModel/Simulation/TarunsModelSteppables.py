from cc3d.core.PySteppables import *
import numpy as np

min_to_mcs = 60.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 15.0

virus_infection_feedback = 1

model_string = '''
// Equations
J1: D -> E; dE*D;
J2: E -> D; dE*E;
J3: E -> Ev; bE*V*E;
J4: Ev -> E; aE*Ev;
J5: Ev -> D; dE*Ev;
J6: Ev -> D; kE*g*Ev*Tc;
J7: -> V; pV*Ev;
J8: V ->; cV*V;
J9: -> Da; bD*V*(D0-Da); //D0 and Da are APC in lung
J10: Da ->; dD*Da;
J11: -> Dm; kD*Da; // Dm are apc in lymph
J12: Dm ->; dDm*Dm;
J13: -> Tc; dc*Tc0;
J14: Tc ->; dc*Tc;
J15: Dm -> Tc; 0.0 * pT1*Dm*Tc/(Dm+pT2) + Dm;
J16: Tc ->; 0.0 * dT1*Tc*Ev/(Ev+dT2);
J17: Th2 -> Th1; (s1*Th1)/(1+Th2)^2 +Th2;
J18: Dm -> Th1; 0.0 * p1*((Da+Dm)*Th1^2)/(1+Th2)^2 + Dm;
J19: Th1 ->; d1*((Da+Dm)*Th1^3)/(1+Th2);
J20: Th1 ->; m*Th1;
J21: -> Th2; s2*Th2/(1+Th2);
J22: Th1 -> Th2; p2*((ro+Th1)*(Da+Dm)*Th1^2)/((1+Th2)*(1+Th2+Th1)) + Th1;
J23: Th2 ->; m*Th2;

// Parameters
aE=5.0*10^-2; 
bE=2*10^-5; 
dE=10^-3;
kE=1.19*10^-2 / 150.0; 
g=0.15 * 150.0; 
tC=0.5;
pV=1.9; 
cV=1;
D0=10^3; 
bD=10^-2; 
dD=2.9;
kD=0.5; 
tD=10; 
dDm=0.5;
dc=1.5*10^-3; 
pT1=2.7; 
pT2=600; 
dT1=2; 
dT2=1;
s1=1.1; 
p1=2; 
d1=0.1; 
m=0.25;
s2=0.1; 
p2=0.01; 
ro=0.5;

// Initial Conditions
E0= 5.0*10^5; 
Tc0 = 1*10^3;
E = E0;
V = 10;
Th1=100;
Th2=100;
'''

lymph_node_string = '''
J11: -> Dm; kD*Da; // Dm are apc in lymph
J12: Dm ->; dDm*Dm;
J13: -> Tc; dc*Tc0;
J14: Tc ->; dc*Tc;
J15: Dm -> Tc; pT1*Dm*Tc/(Dm+pT2) + Dm;
J16: Tc ->; dT1*Tc*Ev/(Ev+dT2);
J17: Th2 -> Th1; (s1*Th1)/(1+Th2)^2 +Th2;
J18: Dm -> Th1; ( p1*((Da+Dm)*Th1^2)/(1+Th2)^2 + Dm);
J19: Th1 ->; d1*((Da+Dm)*Th1^3)/(1+Th2);
J20: Th1 ->; m*Th1;
J21: -> Th2; s2*Th2/(1+Th2);
J22: Th1 -> Th2; p2*((ro+Th1)*(Da+Dm)*Th1^2)/((1+Th2)*(1+Th2+Th1)) + Th1;
J23: Th2 ->; m*Th2;

// Parameters
aE=5.0*10^-2; 
bE=2*10^-5; 
dE=10^-3;
kE=1.19*10^-2 / 150.0; 
g=0.15 * 150.0; 
tC=0.5;
pV=1.9; 
cV=1;
D0=10^3; 
bD=10^-2; 
dD=2.9;
kD=0.5; 
tD=10; 
dDm=0.5;
dc=1.5*10^-3; 
pT1=2.7; 
pT2=600; 
dT1=2; 
dT2=1;
s1=1.1; 
p1=2; 
d1=0.1; 
m=0.25;
s2=0.1; 
p2=0.01; 
ro=0.5;
// Initial Conditions
Tc0 = 1*10^3;
Th1=100;
Th2=100;
// inputs
Da=0.0; 
Ev = 0.0;

// outputs

// Tc

'''


class TarunsModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Updating max simulation steps using scaling factor to simulate 10 days
        self.add_free_floating_antimony(model_string=model_string, model_name='FullModel', step_size=days_to_mcs)
        self.add_free_floating_antimony(model_string=lymph_node_string, model_name='LymphModel', step_size=days_to_mcs)
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs
        self.initial_uninfected = len(self.cell_list_by_type(self.E))
        self.get_xml_element('virus_decay').cdata = self.sbml.FullModel['cV'] * days_to_mcs
        self.scalar_virus = self.sbml.FullModel['V']
        self.shared_steppable_vars['CC3D_lymph_APC_count'] = 0

        numberAPC = 0
        while numberAPC < round(self.sbml.FullModel['D0'] / self.sbml.FullModel['E0'] * self.initial_uninfected):
            x = np.random.randint(10, self.dim.x - 10)
            y = np.random.randint(10, self.dim.y - 10)
            if not self.cell_field[x, y, 0]:
                cell = self.new_cell(self.APC)
                self.cell_field[x:x + 3, y:y + 3, 0] = cell
                cell.targetVolume = cell.volume
                cell.lambdaVolume = cell.volume
                cell.dict['Activation_State'] = False
                numberAPC += 1

    def J1_DtoE(self, cell):
        ## Transition from D to E
        # J1: D -> E; dE*D;
        if self.sbml.FullModel['dE'] * days_to_mcs > np.random.random():
            cell.type = self.E

    def J2_EtoD(self, cell):
        ## Transition from E to D
        # J2: E -> D; dE * E;
        if self.sbml.FullModel['dE'] * days_to_mcs > np.random.random():
            cell.type = self.D

    def J3_EtoEv(self, cell, secretor):
        ## Transition from E to Ev
        # J3: E -> Ev; bE*V*E;
        if virus_infection_feedback == 1:
            bE = self.sbml.FullModel['bE'] * days_to_mcs
            V = self.sbml.FullModel['V']
        elif virus_infection_feedback == 2:
            bE = self.sbml.FullModel['bE'] * days_to_mcs
            V = self.scalar_virus
        elif virus_infection_feedback == 3:
            bE = self.sbml.FullModel['bE'] * days_to_mcs * self.initial_uninfected
            V = secretor.amountSeenByCell(cell)
        else:  # in case of things breaking have a default
            bE = self.sbml.FullModel['bE'] * days_to_mcs
            V = self.sbml.FullModel['V']
        p_EtoEv = bE * V
        if p_EtoEv > np.random.random():
            cell.type = self.EV

    def J4_EvtoE(self, cell):
        ## Transition from Ev to E
        # J4: Ev -> E; aE*Ev;
        if self.sbml.FullModel['aE'] * days_to_mcs > np.random.random():
            cell.type = self.E

    def J5_EvtoD(self, cell):
        ## Transition from Ev to D
        # J5: Ev -> D; dE*Ev;
        dE = self.sbml.FullModel['dE'] * days_to_mcs
        p_EvtoD = dE
        if p_EvtoD > np.random.random():
            cell.type = self.D

    def J6_EvtoD(self, cell, Tc):
        ## Transition from Ev to D
        # J6: Ev -> D; kE * g * Ev * Tc;
        kE = self.sbml.FullModel['kE'] * days_to_mcs
        g = self.sbml.FullModel['g']
        p_EvtoD = kE * g * Tc
        if p_EvtoD > np.random.random():
            cell.type = self.D

    def step(self, mcs):
        for cell in self.cell_list_by_type(self.D):
            self.J1_DtoE(cell)

        for cell in self.cell_list_by_type(self.E):
            self.J2_EtoD(cell)

        secretor = self.get_field_secretor("Virus")
        for cell in self.cell_list_by_type(self.E):
            self.J3_EtoEv(cell, secretor)

        for cell in self.cell_list_by_type(self.EV):
            self.J4_EvtoE(cell)

        for cell in self.cell_list_by_type(self.EV):
            self.J5_EvtoD(cell)

        # TODO: this rescaling needs to be revised. Is matching close enough
        # Tc = len(self.cell_list_by_type(self.TCELL))
        Tc = self.sbml.FullModel['Tc']
        for cell in self.cell_list_by_type(self.EV):
            self.J6_EvtoD(cell, Tc)

        virus_production = 0.0
        for cell in self.cell_list_by_type(self.EV):
            ## Virus Production
            # J7: -> V; pV*Ev;
            pV = self.sbml.FullModel['pV'] * days_to_mcs / self.initial_uninfected * self.sbml.FullModel['E0']
            release = secretor.secreteInsideCellTotalCount(cell, pV / cell.volume)
            virus_production += abs(release.tot_amount)

        ## Virus Decay
        # J8: V ->; cV*V;
        cV = self.sbml.FullModel['cV'] * days_to_mcs
        virus_decay = cV * self.scalar_virus
        self.scalar_virus += virus_production - virus_decay
        self.shared_steppable_vars['scalar_virus'] = self.scalar_virus
        activated_APC_count = 0

        for cell in self.cell_list_by_type(self.APC):
            if not cell.dict['Activation_State']:
                ## Infection and Activation of APC
                # J9: -> Da; bD*V*(D0-Da);
                bD = (self.sbml.FullModel['bD'] * days_to_mcs) * (self.sbml.FullModel['D0'] * self.initial_uninfected /
                                                                  self.sbml.FullModel['E0'])
                # V should be local instead of the total virus
                V = secretor.amountSeenByCell(cell) * self.initial_uninfected
                # V = self.sbml.FullModel['V'] / self.sbml.FullModel['E0'] * self.initial_uninfected
                p_DtoAPC = bD * V
                if p_DtoAPC > np.random.random():
                    cell.dict['Activation_State'] = True
                    activated_APC_count += 1

            elif cell.dict['Activation_State']:
                ## Deactivation of APC
                # J10: Da -> ; dD*Da;
                activated_APC_count += 1
                dD = self.sbml.FullModel['dD'] * days_to_mcs
                p_APCtoD = dD
                if p_APCtoD > np.random.random():
                    cell.dict['Activation_State'] = False
                    activated_APC_count -= 1

        ## Tcell seeding
        # J13: -> Tc; dc*g*Tc0;
        dc = self.sbml.FullModel['dc'] * days_to_mcs
        g = self.sbml.FullModel['g']
        Tc0 = self.sbml.FullModel['Tc0'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        cell = False
        if dc * g * Tc0 > np.random.random():
            while not cell:
                x = np.random.randint(10, self.dim.x - 10)
                y = np.random.randint(10, self.dim.y - 10)
                if not self.cell_field[x, y, 0]:
                    cell = self.new_cell(self.TCELL)
                    self.cell_field[x:x + 3, y:y + 3, 0] = cell
                    cell.targetVolume = cell.volume
                    cell.lambdaVolume = cell.volume
                    cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
                    cd.setLambda(0)
                    cd.assignChemotactTowardsVectorTypes([self.MEDIUM])

        ## Clearance of Tcells
        # J14: Tc ->; dc*Tc;
        for cell in self.cell_list_by_type(self.TCELL):
            dc = self.sbml.FullModel['dc'] * days_to_mcs
            if dc > np.random.random():
                cell.targetVolume = 0.0

        ## Tcell seeding
        # J15a: Dm -> Tc; Dm ;
        Dm = self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        print(Dm)
        if Dm > np.random.random():
            cells_to_seed = max(1, round(Dm))
            print('cells_to_seed', cells_to_seed)
            for i in range(cells_to_seed):
                cell = False
                while not cell:
                    x = np.random.randint(10, self.dim.x - 10)
                    y = np.random.randint(10, self.dim.y - 10)
                    if not self.cell_field[x, y, 0]:
                        cell = self.new_cell(self.TCELL)
                        self.cell_field[x:x + 3, y:y + 3, 0] = cell
                        cell.targetVolume = cell.volume
                        cell.lambdaVolume = cell.volume
                        cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
                        cd.setLambda(0)
                        cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
                        print('Cell')

        ## Step SBML forward
        self.timestep_sbml()


class ChemotaxisSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        for cell in self.cell_list_by_type(self.APC, self.TCELL):
            cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
            cd.setLambda(0)
            cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
        self.secretor = self.get_field_secretor("Virus")

    def step(self, mcs):
        lambda_chemotaxis = 250.0
        for cell in self.cell_list_by_type(self.APC, self.TCELL):
            cd = self.chemotaxisPlugin.getChemotaxisData(cell, "Virus")
            cd.setLambda(0)
            concentration = self.secretor.amountSeenByCell(cell)
            if concentration:
                cd.setLambda(lambda_chemotaxis / (1 + concentration))


class PlotsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.seed = 0
        self.initial_uninfected = len(self.cell_list_by_type(self.E))

        self.plot_win = self.add_new_plot_window(title='Epithelial Cells',
                                                 x_axis_title='Time (days)',
                                                 y_axis_title='Fraction of Cells', x_scale_type='linear',
                                                 y_scale_type='linear',
                                                 grid=False)

        self.plot_win.add_plot("ODEE", style='Dots', color='blue', size=5)
        self.plot_win.add_plot("ODEEv", style='Dots', color='yellow', size=5)
        self.plot_win.add_plot("ODED", style='Dots', color='red', size=5)
        self.plot_win.add_plot("CC3DE", style='Lines', color='blue', size=5)
        self.plot_win.add_plot("CC3DEv", style='Lines', color='yellow', size=5)
        self.plot_win.add_plot("CC3DD", style='Lines', color='red', size=5)

        self.plot_win2 = self.add_new_plot_window(title='Virus',
                                                  x_axis_title='Time (days)',
                                                  y_axis_title='Virus', x_scale_type='linear', y_scale_type='log',
                                                  grid=False)

        self.plot_win2.add_plot("ODEV", style='Dots', color='red', size=5)
        self.plot_win2.add_plot("CC3DV", style='Lines', color='red', size=5)

        self.plot_win3 = self.add_new_plot_window(title='Epithelial Tcells',
                                                  x_axis_title='Time (days)',
                                                  y_axis_title='Number of Cells', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)

        self.plot_win3.add_plot("ODETc", style='Dots', color='red', size=5)
        self.plot_win3.add_plot("CC3DTc", style='Lines', color='red', size=5)

        self.plot_win4 = self.add_new_plot_window(title='Epithelial APC',
                                                  x_axis_title='Time (days)',
                                                  y_axis_title='Number of Cells', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)

        self.plot_win4.add_plot("ODEAPC", style='Dots', color='red', size=5)
        self.plot_win4.add_plot("CC3DAPC0", style='Lines', color='orange', size=5)
        self.plot_win4.add_plot("CC3DAPC", style='Lines', color='red', size=5)

        self.plot_win5 = self.add_new_plot_window(title='Lymph node APC',
                                                  x_axis_title='Time (days)',
                                                  y_axis_title='Number of Cells', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)

        self.plot_win5.add_plot("ODEAPC", style='Dots', color='red', size=5)
        # self.plot_win5.add_plot("CC3DAPC0", style='Lines', color='orange', size=5)
        self.plot_win5.add_plot("CC3DAPC", style='Lines', color='red', size=5)

    def step(self, mcs):
        self.plot_win.add_data_point("ODEE", mcs * days_to_mcs, self.sbml.FullModel['E'] / self.sbml.FullModel['E0'])
        self.plot_win.add_data_point("ODEEv", mcs * days_to_mcs, self.sbml.FullModel['Ev'] / self.sbml.FullModel['E0'])
        self.plot_win.add_data_point("ODED", mcs * days_to_mcs, self.sbml.FullModel['D'] / self.sbml.FullModel['E0'])

        self.plot_win2.add_data_point("ODEV", mcs * days_to_mcs, self.sbml.FullModel['V'])
        self.plot_win3.add_data_point("ODETc", mcs * days_to_mcs, self.sbml.FullModel['g'] *
                                      self.sbml.FullModel['Tc'] / self.sbml.FullModel['E0'] * self.initial_uninfected)
        self.plot_win4.add_data_point("ODEAPC", mcs * days_to_mcs,
                                      self.sbml.FullModel['Da'] / self.sbml.FullModel['E0'])

        self.plot_win5.add_data_point("ODEAPC", mcs * days_to_mcs,
                                      self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'])

        secretor = self.get_field_secretor("Virus")
        self.field_virus = 0.0
        for cell in self.cell_list:
            self.field_virus += secretor.amountSeenByCell(cell)

        self.num_activeAPC = 0.0
        for cell in self.cell_list_by_type(self.APC):
            if cell.dict['Activation_State']:
                self.num_activeAPC += 1

        self.plot_win.add_data_point("CC3DE", mcs * days_to_mcs,
                                     len(self.cell_list_by_type(self.E)) / self.initial_uninfected)
        self.plot_win.add_data_point("CC3DEv", mcs * days_to_mcs,
                                     len(self.cell_list_by_type(self.EV)) / self.initial_uninfected)
        self.plot_win.add_data_point("CC3DD", mcs * days_to_mcs,
                                     len(self.cell_list_by_type(self.D)) / self.initial_uninfected)
        self.plot_win2.add_data_point("CC3DV", mcs * days_to_mcs,
                                      self.field_virus)
        self.plot_win4.add_data_point("CC3DAPC0", mcs * days_to_mcs,
                                      len(self.cell_list_by_type(self.APC)) / self.initial_uninfected)
        self.plot_win4.add_data_point("CC3DAPC", mcs * days_to_mcs,
                                      self.num_activeAPC / self.initial_uninfected)
        self.plot_win3.add_data_point("CC3DTc", mcs * days_to_mcs,
                                      len(self.cell_list_by_type(self.TCELL)))
        self.plot_win5.add_data_point("CC3DAPC", mcs * days_to_mcs,
                                      self.shared_steppable_vars['CC3D_lymph_APC_count'] / self.initial_uninfected)
