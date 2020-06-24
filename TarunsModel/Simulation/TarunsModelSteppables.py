from cc3d.core.PySteppables import *
import numpy as np

min_to_mcs = 60.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 10.0

virus_infection_feedback = 1

model_string = '''
#reactions
J1: D -> E; dE*D;
J2: E -> D; dE*E;
J3: E -> Ev; bE*V*E;
J4: Ev -> E; aE*Ev;
J5: Ev -> D; dE*Ev;
J6: Ev -> D; kE*g*Ev*Tc;
J7: -> V; pV*Ev;
J8: V ->; cV*V;
J9: -> Da; bD*V*(D0-Da);
J10: Da ->; dD*Da;
J11: -> Dm; kD*Da;
J12: Dm ->; dDm*Dm;
J13: -> Tc; dc*Tc0;
J14: Tc ->; dc*Tc;
J15: Dm -> Tc; pT1*Dm*Tc/(Dm+pT2) + Dm;
J16: Tc ->; dT1*Tc*Ev/(Ev+dT2);
J17: Th2 -> Th1; (s1*Th1)/(1+Th2)^2 +Th2;
J18: Dm -> Th1; p1*((Da+Dm)*Th1^2)/(1+Th2)^2 + Dm;
J19: Th1 ->; d1*((Da+Dm)*Th1^3)/(1+Th2);
J20: Th1 ->; m*Th1;
J21: -> Th2; s2*Th2/(1+Th2);
J22: Th1 -> Th2; p2*((ro+Th1)*(Da+Dm)*Th1^2)/((1+Th2)*(1+Th2+Th1)) + Th1;
J23: Th2 ->; m*Th2;

//parameters
E0=5.0*10^5; 
aE=5.0*10^-2; 
bE=2*10^-5; 
dE=10^-3;
kE=1.19*10^-2; 
g=0.15; 
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
Tc0=1*10^3; 
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
E = E0;
V = 10;
Th1=100;
Th2=100;
'''

class TarunsModelSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        # Uptading max simulation steps using scaling factor to simulate 10 days
        self.add_free_floating_antimony(model_string=model_string, model_name='FullModel', step_size=days_to_mcs)
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs
        self.initial_uninfected = len(self.cell_list_by_type(self.E))
        self.get_xml_element('virus_decay').cdata = self.sbml.FullModel['cV'] * days_to_mcs
        self.scalar_virus = self.sbml.FullModel['V']

    def step(self,mcs):
        secretor = self.get_field_secretor("Virus")
        for cell in self.cell_list_by_type(self.E):
            ## Transition from E to Ev
            # J3: E -> Ev; bE*V*E;
            if virus_infection_feedback == 1:
                bE = self.sbml.FullModel['bE'] * days_to_mcs
                V = self.sbml.FullModel['V']
            if virus_infection_feedback == 2:
                bE = self.sbml.FullModel['bE'] * days_to_mcs
                V = self.scalar_virus
            if virus_infection_feedback == 3:
                bE = self.sbml.FullModel['bE'] * days_to_mcs * self.initial_uninfected
                V = secretor.amountSeenByCell(cell)
            p_EtoEv = bE * V
            if p_EtoEv > np.random.random():
                cell.type = self.EV

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

        ## Step SBML forward
        self.timestep_sbml()

class PlotsSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):
        self.initial_uninfected = len(self.cell_list_by_type(self.E))

        self.plot_win = self.add_new_plot_window(title='Fraction of Cells',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Fraction of Uninfected Cells', x_scale_type='linear',
                                                 y_scale_type='linear',
                                                 grid=False)

        self.plot_win.add_plot("ODEE", style='Dots', color='blue', size=5)
        self.plot_win.add_plot("CC3DE", style='Lines', color='blue', size=5)

        self.plot_win2 = self.add_new_plot_window(title='Virus',
                                                  x_axis_title='MonteCarlo Step (MCS)',
                                                  y_axis_title='Virus', x_scale_type='linear', y_scale_type='log',
                                                  grid=False)

        self.plot_win2.add_plot("ODEV", style='Dots', color='red', size=5)
        self.plot_win2.add_plot("CC3DVS", style='Lines', color='red', size=5)
        self.plot_win2.add_plot("CC3DVF", style='Lines', color='orange', size=5)

    def step(self, mcs):
        self.scalar_virus = self.shared_steppable_vars['scalar_virus']

        self.plot_win.add_data_point("ODEE", mcs * days_to_mcs,self.sbml.FullModel['E']/self.sbml.FullModel['E0'])
        self.plot_win2.add_data_point("ODEV", mcs * days_to_mcs,self.sbml.FullModel['V'])

        self.plot_win.add_data_point("CC3DE", mcs * days_to_mcs,
                                     len(self.cell_list_by_type(self.E))/self.initial_uninfected)
        self.field_virus = 0.0
        secretor = self.get_field_secretor("Virus")
        self.field_virus = 0.0
        for cell in self.cell_list:
            self.field_virus += secretor.amountSeenByCell(cell)
        self.plot_win2.add_data_point("CC3DVS", mcs * days_to_mcs, self.scalar_virus)
        self.plot_win2.add_data_point("CC3DVF", mcs * days_to_mcs, self.field_virus)