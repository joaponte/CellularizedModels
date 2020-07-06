from cc3d.core.PySteppables import *
import numpy as np

min_to_mcs = 60.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 30.0

virus_infection_feedback = 1
use_LymphModel_outputs = False

contact_cytotoxicity = True

plot_Epithelial_cells = True
plot_Virus = True
plot_Tcells = True
plot_APC = True
plot_Dm = True

plot_Contact_Cytotoxicity = True

# TODO: add graphs for all state variables

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

J9: -> Da; bD*V*(D0-Da);
J10: Da ->; dD*Da;

// bellow is lymph

J11: -> Dm; kD*Da;
J12: Dm ->; dDm*Dm;
J13: -> Tc; dC*Tc0;
J14: Tc ->; dC*Tc;
//# J15: Dm -> Tc; rT1*Dm*Tc/(Dm+pT2) + Dm; 
J15:  -> Tc; rT1*Dm*Tc/(Dm+pT2); 
J16: Tc ->; dT1*Tc*Ev/(Ev+dT2);

J17: -> Th1; sTh1*Th1/(1+Th2)^2;
//J18: Dm -> Th1; pTh1*Dm*(Th1^2)/(1+Th2)^2 + Dm;
J18:  -> Th1; pTh1*Dm*(Th1^2)/(1+Th2)^2;
J19: Th1 ->; dTh1*Dm*(Th1^3)/(1+Th2);
J20: Th1 ->; mTh*Th1;
J21: -> Th2; sTh2*Th2/(1+Th2);
# J22: Dm -> Th2; pTh2*(r+Th1)*Dm*(Th2^2)/((1+Th2)*(1+Th1+Th2)) + Dm
J22:  -> Th2; pTh2*(r+Th1)*Dm*(Th2^2)/((1+Th2)*(1+Th1+Th2)) 
J23: Th2 ->; mTh*Th2;

// new eqs

J24: -> B; dB*B0;
J25: B ->; dB*B;
//J26: Dm + Th2 -> B; rB1*B*(Dm+h*Th2)/(Dm+h*Th2+rB2);
J26:  -> B; rB1*B*(Dm+h*Th2)/(Dm+h*Th2+rB2);

J27: B -> Pss; pS*B;
J28: B -> Psn; pS*B;
J29: B -> Pls; pL*B*Th2;
J30: B -> Pln; pL*B*Th2;
J31: Pss ->; dS*Pss;
J32: Psn ->; dS*Psn;
J33: Pls ->; dL*Pls;
J34: Pln ->; dL*Pln;
J35: Pss -> Pls; d*(1-v)*Pss;
J36: Psn -> Pln; d*(1-v)*Psn;
J37:  -> Pss; b*v*Pss;
J37a:  -> Pls; b*v*Pls;
J38:  -> Psn; b*v*Psn;
J38a:  -> Pln; b*v*Pln;
J39: Pls -> Pss; d*(1-v)*Pls;
J40: Pln -> Psn; d*(1-v)*Pln; 

J41:  -> sIgM; pAS*Pss;
J42:  -> nIgM; pAS*Psn;
J43:  -> sIgG; pAS*Pls;
J44:  -> nIgG; pAS*Pln;
J45: sIgM ->; dM*sIgM;
J46: sIgG ->; dG*sIgG;
J47: nIgM ->; dM*nIgM;
J48: nIgG ->; dG*nIgG; 
// feed back to tissue
J49: Ev -> D; eE*Ev*nIgM; // reevaluate these set (49-52) due to anti-bodies being consumed 
J50: Ev -> D; eE*Ev*nIgG;
J51: V  ->; eV*V*sIgM;
//J51: V + sIgM ->; eV*V*sIgM;
J52: V  ->; eV*V*sIgG;
//J52: V + sIgG ->; eV*V*sIgG;


// Parameters
dE=10^-3;
E0=5*10^5;
bE=3*10^-6;
aE=5.0*10^-2;
V0=10;
pV=19;
cV=1;
kE=1.19*10^-3 / 9; // 900 rescaling for non-0 T cell pop in tissue
g=0.15 * 9;
tC=0.5;// not included in the model
eE=0.05;
eV=16;
D0=10^3;
bD=10^-2;
dD=2.9;
kD=0.5;
tD=1;// not included in the model
dDm=0.5;
dC=10.1*10^-3;
Tc0=5*10^2;
rT1=1.3;
rT2=1;
dT1=5.0;
dT2=2000;
sTh1=2.5;
pTh1=4;
dTh1=0.03;
mTh=0.25;
sTh2=0.001;
pTh2=0.0012;
r=0.2;
dB=0.0009;
B0=1*10^3;
rB1=100;
rB2=2*10^5;
h=100;
pS=10^-1;
pL=8*10^-3;
dS=0.002;
dL=0.02;
b=2.4*10^-4;
d=2.4*10^-2;
pAS=0.2;
pAL=0.3;
dG=0.04;
dM=0.2;
pT2=600;
v = 0.5
// Initial Conditions
E = E0;
V = 10;
Th1=10;
Th2=10;
'''

lymph_node_string = '''
J11: -> Dm; 0.0*kD*Da; // MUST STAY SHUT OFF, DM IS AN INPUT Dm are apc in lymph
J12: Dm ->; dDm*Dm;
J13: -> Tc; dC*Tc0;
J14: Tc ->; dC*Tc;
// J15: Dm -> Tc; (rT1*Dm*Tc/(Dm+pT2) + Dm);
J15:  -> Tc; rT1*Dm*Tc/(Dm+pT2); 
J16: Tc ->; dT1*Tc*Ev/(Ev+dT2);
J17: -> Th1; sTh1*Th1/(1+Th2)^2;
//J18: Dm -> Th1; pTh1*Dm*(Th1^2)/(1+Th2)^2 + Dm;
J18:  -> Th1; pTh1*Dm*(Th1^2)/(1+Th2)^2;
J19: Th1 ->; dTh1*Dm*(Th1^3)/(1+Th2);
J20: Th1 ->; mTh*Th1;
J21: -> Th2; sTh2*Th2/(1+Th2);
//J22: Dm -> Th2; pTh2*(r+Th1)*Dm*(Th2^2)/((1+Th2)*(1+Th1+Th2)) + Dm
J22:  -> Th2; pTh2*(r+Th1)*Dm*(Th2^2)/((1+Th2)*(1+Th1+Th2)) 
J23: Th2 ->; mTh*Th2;
// new eqs

J24: -> B; dB*B0;
J25: B ->; dB*B;
//J26: Dm + Th2 -> B; rB1*B*(Dm+h*Th2)/(Dm+h*Th2+rB2);
J26:  -> B; rB1*B*(Dm+h*Th2)/(Dm+h*Th2+rB2);

J27: B -> Pss; pS*B;
J28: B -> Psn; pS*B;
J29: B -> Pls; pL*B*Th2;
J30: B -> Pln; pL*B*Th2;
J31: Pss ->; dS*Pss;
J32: Psn ->; dS*Psn;
J33: Pls ->; dL*Pls;
J34: Pln ->; dL*Pln;
J35: Pss -> Pls; d*(1-v)*Pss;
J36: Psn -> Pln; d*(1-v)*Psn;
J37: Pss -> Pss; b*v*Pss;
J38: Psn -> Psn; b*v*Psn;
J39: Pls -> Pss; d*(1-v)*Pls;
J40: Pln -> Psn; d*(1-v)*Pln; 

J41:  -> sIgM; pAS*Pss;
J42:  -> nIgM; pAS*Psn;
J43:  -> sIgG; pAS*Pls;
J44:  -> nIgG; pAS*Pln;
J45: sIgM ->; dM*sIgM;
J46: sIgG ->; dG*sIgG;
J47: nIgM ->; dM*nIgM;
J48: nIgG ->; dG*nIgG; 

// feed back to tissue
J49: Ev -> D; eE*Ev*nIgM;
J50: Ev -> D; eE*Ev*nIgG;
J51: V ->; eV*V*sIgM;
J52: V ->; eV*V*sIgG;

// Parameters
dE=10^-3;
E0=5*10^5;
bE=3*10^-6;
aE=5.0*10^-2;
V0=10;
pV=19;
cV=1;
kE=1.19*10^-3 / 900; // 900 rescaling for non-0 T cell pop in tissue
g=0.15 * 900;
tC=0.5;// not included in the model
eE=0.05;
eV=16;
D0=10^3;
bD=10^-2;
dD=2.9;
kD=0.5;
tD=1;// not included in the model
dDm=0.5;
dC=10.1*10^-3;
Tc0=5*10^2;
rT1=1.3;
rT2=1;
dT1=5.0;
dT2=2000;
sTh1=2.5;
pTh1=4;
dTh1=0.03;
mTh=0.25;
sTh2=0.001;
pTh2=0.0012;
r=0.2;
dB=0.0009;
B0=1*10^3;
rB1=100;
rB2=2*10^5;
h=100;
pS=10^-1;
pL=8*10^-3;
dS=0.002;
dL=0.02;
b=2.4*10^-4;
d=2.4*10^-2;
pAS=0.2;
pAL=0.3;
dG=0.04;
dM=0.2;
pT2=600;
v = 0.5;

// Initial Conditions
E = E0;
V = 10;
Th1=10;
Th2=10;
// inputs
Da=0.0; // kD*Da
Ev = 0.0;
// ck1

// outputs

// Tc
// ck1
// ck2
'''


class TarunsModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def seed_epithelial_sheet(self):
        for x in range(0, self.dim.x, int(3)):
            for y in range(0, self.dim.y, int(3)):
                cell = self.new_cell(self.E)
                self.cellField[x:x + int(3), y:y + int(3), 0] = cell

    def init_variables(self):
        # Updating max simulation steps using scaling factor to simulate 10 days
        self.add_free_floating_antimony(model_string=model_string, model_name='FullModel', step_size=days_to_mcs)
        self.add_free_floating_antimony(model_string=lymph_node_string, model_name='LymphModel', step_size=days_to_mcs)
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs
        self.initial_uninfected = len(self.cell_list_by_type(self.E))
        self.get_xml_element('virus_decay').cdata = self.sbml.FullModel['cV'] * days_to_mcs
        self.scalar_virus = self.sbml.FullModel['V']
        self.shared_steppable_vars['CC3D_lymph_APC_count'] = 0.0
        self.shared_steppable_vars['Active_APCs'] = 0.0
        self.shared_steppable_vars['ODE_Killing'] = 0.0
        self.shared_steppable_vars['Contact_Killing'] = 0.0
        self.virus_secretor = self.get_field_secretor("Virus")

    def seed_APC_cells(self):
        numberAPC = 0.0
        while numberAPC < round(self.sbml.FullModel['D0'] / self.sbml.FullModel['E0'] * self.initial_uninfected):
            x = np.random.randint(10, self.dim.x - 10)
            y = np.random.randint(10, self.dim.y - 10)
            if not self.cell_field[x, y, 1]:
                cell = self.new_cell(self.APC)
                self.cell_field[x:x + 3, y:y + 3, 1] = cell
                cell.targetVolume = cell.volume
                cell.lambdaVolume = cell.volume
                cell.dict['Activation_State'] = 0  # in tissue
                numberAPC += 1

    def start(self):
        self.seed_epithelial_sheet()

        self.init_variables()

        self.seed_APC_cells()

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

    def J3_EtoEv(self, cell):
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
            V = self.virus_secretor.amountSeenByCell(cell)
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
        p_EvtoD = kE * Tc
        if p_EvtoD > np.random.random():
            if not contact_cytotoxicity:
                cell.type = self.D
            self.shared_steppable_vars['ODE_Killing'] += 1

    def J7_virus_production(self):
        virus_production = 0.0
        for cell in self.cell_list_by_type(self.EV):
            ## Virus Production
            # J7: -> V; pV*Ev;
            pV = self.sbml.FullModel['pV'] * days_to_mcs / self.initial_uninfected * self.sbml.FullModel['E0']
            release = self.virus_secretor.secreteInsideCellTotalCount(cell, pV / cell.volume)
            virus_production += abs(release.tot_amount)
        return virus_production

    def J8_virus_decay(self):
        # no longer used, now update_scalar_virus does all the work
        ## Virus Decay
        # J8: V ->; cV*V;
        cV = self.sbml.FullModel['cV'] * days_to_mcs
        virus_decay = cV * self.scalar_virus



    def J9_J10_APC_activation_deactivation(self):
        tissue_apc = 0
        lymph_apc = 0
        node_apc = 0
        just_moved_in_to_node = 0
        # activated_APC_count = 0
        for cell in self.cell_list_by_type(self.APC):
            # if cell.dict == 0 -> in tissue; keep logic (J9: -> Da; bD*V*(D0-Da);) but replace cell.dict[
            # 'Activation_State'] = True by cell.dict['Activation_State'] = 1, this cell is now in transit:
            # tissue_apc -= 1; lymph_apc+=1
            ########################
            # add this elif elif cell.dict['Activation_State'] = 1; go from the lymph to the
            # node by p_transition = 1/delay [rescaled] if p_transition > dice_roll:  cell.dict['Activation_State'] =
            # 2 this cell is now in the lymph; lymph_apc-=1 node_apc+=1 ## replace the elif cell.dict[
            # 'Activation_State']: elif cell.dict['Activation_State'] == 2; keep logic replace cell.dict[
            # 'Activation_State'] = False by cell.dict['Activation_State'] = 0 ## then number of apc in node is
            # cell.dict['Activation_State'] == 2 * self.sbml.LymphModel['kD'] / delay
            if cell.dict['Activation_State'] == 0:  # in tissue
                ## Infection and Activation of APC
                # J9: -> Da; bD*V*(D0-Da);
                tissue_apc += 1
                bD = (self.sbml.FullModel['bD'] * days_to_mcs) * (self.sbml.FullModel['D0'] * self.initial_uninfected /
                                                                  self.sbml.FullModel['E0'])
                # V should be local instead of the total virus
                V = self.virus_secretor.amountSeenByCell(cell) * self.initial_uninfected
                # V = self.sbml.FullModel['V'] / self.sbml.FullModel['E0'] * self.initial_uninfected

                p_tissue_2_lymph = bD * V
                if p_tissue_2_lymph > np.random.random():
                    cell.dict['Activation_State'] = 1
                    tissue_apc -= 1
                    lymph_apc += 1
            elif cell.dict['Activation_State'] == 1:  # in lymph
                # transport of apc to lymph
                lymph_apc += 1
                p_lymph_2_node = 1 / 10  # PLACEHOLDER
                if p_lymph_2_node > np.random.random():
                    cell.dict['Activation_State'] = 2
                    lymph_apc -= 1
                    node_apc += 1
                    just_moved_in_to_node += 1
            elif cell.dict['Activation_State'] == 2:
                node_apc += 1
                dD = self.sbml.FullModel['dD'] * days_to_mcs
                p_lymph_2_tissue = dD
                if p_lymph_2_tissue > np.random.random():
                    cell.dict['Activation_State'] = 0
                    node_apc -= 1
                    tissue_apc += 1
        self.shared_steppable_vars['Active_APCs'] = lymph_apc + node_apc
        return tissue_apc, lymph_apc, node_apc, just_moved_in_to_node

    def J11_APC_travel_Lymph_Model_input(self, just_moved_in_to_node):
        ## APC "travel" to lymph node
        # we'll implement it as a signal that is proportional to the activated apcs
        # J11: -> Dm; kD * Da; // Dm are apc in lymph
        # print(node_apc)
        self.sbml.LymphModel['Dm'] += max(0, (just_moved_in_to_node * self.sbml.FullModel['D0'] *
                                              self.initial_uninfected / self.sbml.FullModel['E0']))

    def J13_Tcell_stable_population_seeding(self):
        ## Tcell seeding
        # J13: -> Tc; dc*g*Tc0;
        # TODO: replace FullModel with LymphModel where appropriate
        dc = self.sbml.FullModel['dC'] * days_to_mcs
        g = self.sbml.FullModel['g']
        Tc0 = self.sbml.FullModel['Tc0'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        cell = False
        if dc * g * Tc0 > np.random.random():
            while not cell:
                x = np.random.randint(0, self.dim.x - 3)
                y = np.random.randint(0, self.dim.y - 3)
                if not self.cell_field[x, y, 1]:
                    cell = self.new_cell(self.TCELL)
                    self.cell_field[x:x + 3, y:y + 3, 1] = cell
                    cell.targetVolume = cell.volume
                    cell.lambdaVolume = cell.volume
                    cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
                    cd.setLambda(0)
                    cd.assignChemotactTowardsVectorTypes([self.MEDIUM])

    def J14_Tcell_clearance(self):
        ## Clearance of Tcells
        # J14: Tc ->; dc*Tc;
        for cell in self.cell_list_by_type(self.TCELL):
            dc = self.sbml.FullModel['dC'] * days_to_mcs
            if dc > np.random.random():
                cell.targetVolume = 0.0

    def J15a_Tcell_inflamatory_seeding(self):
        ## Tcell seeding
        # J15a: Dm -> Tc; Dm ;
        g = self.sbml.FullModel['g']
        # TODO: replace FullModel with LymphModel where appropriate

        Dm = self.sbml.FullModel['Dm'] * days_to_mcs / self.sbml.FullModel['E0'] * self.initial_uninfected
        if Dm * g > np.random.random():
            cells_to_seed = max(1, round(Dm * g))
            for i in range(cells_to_seed):
                cell = False
                while not cell:
                    x = np.random.randint(0, self.dim.x - 3)
                    y = np.random.randint(0, self.dim.y - 3)
                    if not self.cell_field[x, y, 1]:
                        cell = self.new_cell(self.TCELL)
                        self.cell_field[x:x + 3, y:y + 3, 1] = cell
                        cell.targetVolume = cell.volume
                        cell.lambdaVolume = cell.volume
                        cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
                        cd.setLambda(0)
                        cd.assignChemotactTowardsVectorTypes([self.MEDIUM])

    def J15b_Tcell_inflamatory_seeding(self):
        ## Tcell seeding
        # J15b: Dm -> Tc; rT1 * Dm * Tc / (Dm + pT2)
        # TODO: replace FullModel with LymphModel where appropriate

        rT1 = self.sbml.FullModel['rT1'] * days_to_mcs
        Dm = self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        Tc = self.sbml.FullModel['Tc'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        pT2 = self.sbml.FullModel['pT2'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        g = self.sbml.FullModel['g']
        if g * rT1 * Dm * Tc / (Dm + pT2) > np.random.random():
            cells_to_seed = max(1, round(g * rT1 * Dm * Tc / (Dm + pT2)))
            for i in range(cells_to_seed):
                cell = False
                while not cell:
                    x = np.random.randint(0, self.dim.x - 3)
                    y = np.random.randint(0, self.dim.y - 3)
                    if not self.cell_field[x, y, 1]:
                        cell = self.new_cell(self.TCELL)
                        self.cell_field[x:x + 3, y:y + 3, 1] = cell
                        cell.targetVolume = cell.volume
                        cell.lambdaVolume = cell.volume
                        cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
                        cd.setLambda(0)
                        cd.assignChemotactTowardsVectorTypes([self.MEDIUM])

    def J15_Tcell_inflamatory_seeding(self):
        self.J15a_Tcell_inflamatory_seeding()
        self.J15b_Tcell_inflamatory_seeding()

    def J16_Tcell_clearance(self):
        # TODO: NEEDs RESCALING
        # Tcell clearance
        # J16: Tc ->; dT1 * Tc * Ev/(Ev+dT2)
        Ev = len(self.cell_list_by_type(self.EV))
        dT1 = self.sbml.FullModel['dT1'] * days_to_mcs * 3
        dT2 = self.sbml.FullModel['dT1'] * days_to_mcs / self.sbml.FullModel['E0'] * self.initial_uninfected
        for cell in self.cell_list_by_type(self.TCELL):
            if dT1 * Ev / (Ev + dT2) > np.random.random():
                cell.targetVolume = 0.0

    def J51_virus_decay_from_sIgM(self):
        # J51: V ->; eV * V * sIgM;

        if use_LymphModel_outputs:
            sIgM = self.sbml.LymphModel['sIgM']
        else:
            sIgM = self.sbml.FullModel['sIgM']
        eV = self.sbml.FullModel['eV']
        gamma_sIgM = eV * sIgM
        return gamma_sIgM

    def J52_virus_decay_from_sIgG(self):
        # J52: V ->; eV * V * sIgG;
        if use_LymphModel_outputs:
            sIgG = self.sbml.LymphModel['sIgG']
        else:
            sIgG = self.sbml.FullModel['sIgG']
        eV = self.sbml.FullModel['eV']
        gamma_sIgG = eV * sIgG
        return gamma_sIgG

    def J51_J52_update_virus_decay(self):

        gamma_sIgM = self.J51_virus_decay_from_sIgM()
        gamma_sIgG = self.J52_virus_decay_from_sIgG()
        # print((self.sbml.FullModel['cV'] + gamma_sIgM) * days_to_mcs)
        self.get_xml_element('virus_decay').cdata = min(.999,
                                                        (self.sbml.FullModel['cV'] + gamma_sIgM + gamma_sIgG) *
                                                        days_to_mcs)

        return min(.999, (self.sbml.FullModel['cV'] + gamma_sIgM + gamma_sIgG) * days_to_mcs)

    def update_scalar_virus(self, virus_production, current_gamma):
        # does the work that J8 used to do
        virus_decay = current_gamma * self.scalar_virus
        self.scalar_virus += virus_production - virus_decay
        self.shared_steppable_vars['scalar_virus'] = self.scalar_virus

    def lymph_model_input_from_full(self, Ev, Da):
        Ev *= self.sbml.FullModel['E0'] / self.initial_uninfected
        self.sbml.LymphModel['Ev'] = Ev
        self.sbml.LymphModel['Da'] = Da
        self.sbml.LymphModel['V'] = self.scalar_virus

    def contact_killing(self, cell):

        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
            if neighbor:
                if neighbor.type == self.EV:
                    self.shared_steppable_vars['Contact_Killing'] += 1
                    if contact_cytotoxicity:
                        neighbor.type = self.D

    def test_prints(self):
        print('taruns\tcc3d')
        print('Tc', self.sbml.FullModel['Tc'], self.sbml.LymphModel['Tc'])
        print('Th1', self.sbml.FullModel['Th1'], self.sbml.LymphModel['Th1'])
        print('Th2', self.sbml.FullModel['Th2'], self.sbml.LymphModel['Th2'])
        # print('Da', self.sbml.FullModel['Da'], self.sbml.LymphModel['Da'])
        print('Dm', self.sbml.FullModel['Dm'], self.sbml.LymphModel['Dm'])

    def step(self, mcs):
        # self.test_prints()
        self.lymph_model_input_from_full(len(self.cell_list_by_type(self.EV)), self.sbml.FullModel['Da'])
        for cell in self.cell_list_by_type(self.D):
            self.J1_DtoE(cell)

        for cell in self.cell_list_by_type(self.E):
            self.J2_EtoD(cell)

        # self.virus_secretor = self.get_field_secretor("Virus")
        for cell in self.cell_list_by_type(self.E):
            self.J3_EtoEv(cell)

        for cell in self.cell_list_by_type(self.EV):
            self.J4_EvtoE(cell)

        for cell in self.cell_list_by_type(self.EV):
            self.J5_EvtoD(cell)

        g = self.sbml.FullModel['g']
        # Tc = self.sbml.FullModel['Tc'] * g
        Tc = len(self.cell_list_by_type(self.TCELL)) / self.initial_uninfected * self.sbml.FullModel['E0']
        for cell in self.cell_list_by_type(self.EV):
            self.J6_EvtoD(cell, Tc)

        virus_production = self.J7_virus_production()

        current_gamma = self.J51_J52_update_virus_decay()
        self.update_scalar_virus(virus_production, current_gamma)
        # activated_APC_count = 0

        tissue_apc, lymph_apc, node_apc, just_moved_in_to_node = self.J9_J10_APC_activation_deactivation()

        self.J11_APC_travel_Lymph_Model_input(just_moved_in_to_node)

        self.J13_Tcell_stable_population_seeding()

        self.J14_Tcell_clearance()

        self.J15_Tcell_inflamatory_seeding()

        self.J16_Tcell_clearance()

        ## Tcell Contact Killing
        for cell in self.cell_list_by_type(self.TCELL):
            self.contact_killing(cell)

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

        if plot_Epithelial_cells:
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

        if plot_Virus:
            self.plot_win2 = self.add_new_plot_window(title='Virus',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Virus', x_scale_type='linear', y_scale_type='log',
                                                      grid=False)

            self.plot_win2.add_plot("ODEV", style='Dots', color='red', size=5)
            self.plot_win2.add_plot("CC3DV", style='Lines', color='red', size=5)

        if plot_Tcells:
            self.plot_win3 = self.add_new_plot_window(title='Epithelial Tcells',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False)

            self.plot_win3.add_plot("ODETc", style='Dots', color='red', size=5)
            self.plot_win3.add_plot("CC3DTc", style='Lines', color='red', size=5)

        if plot_APC:
            self.plot_win4 = self.add_new_plot_window(title='Epithelial APC',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False)

            self.plot_win4.add_plot("ODEAPC", style='Dots', color='red', size=5)
            self.plot_win4.add_plot("CC3DAPC0", style='Lines', color='orange', size=5)
            self.plot_win4.add_plot("CC3DAPC", style='Lines', color='red', size=5)

        if plot_Dm:
            self.plot_win5 = self.add_new_plot_window(title='Lymph node APC',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False)

            self.plot_win5.add_plot("ODEAPC", style='Dots', color='red', size=5)
            # self.plot_win5.add_plot("CC3DAPC0", style='Lines', color='orange', size=5)

            self.plot_win5.add_plot("CC3DAPC", style='Lines', color='red', size=5)

        if (contact_cytotoxicity == True) and (plot_Contact_Cytotoxicity == True):
            self.plot_win6 = self.add_new_plot_window(title='Contact Cytotoxicity',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=False)

            self.plot_win6.add_plot("ODECC", style='Dots', color='red', size=5)
            self.plot_win6.add_plot("CC3DCC", style='Lines', color='red', size=5)

        self.plot_win7 = self.add_new_plot_window(title='Anti-Bodies',
                                                  x_axis_title='Time (days)',
                                                  y_axis_title='Conc', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)

        self.plot_win7.add_plot("nIgM SBML", style='Dots', color='red', size=5)
        self.plot_win7.add_plot("nIgG SBML", style='Lines', color='blue', size=5)

        self.plot_win7.add_plot("sIgM SBML", style='Lines', color='yellow', size=5)
        self.plot_win7.add_plot("sIgG SBML", style='Lines', color='magenta', size=5)

    def step(self, mcs):

        # self.plot_win7.add_data_point("nIgM SBML", mcs * days_to_mcs, self.sbml.FullModel['nIgM'])
        # self.plot_win7.add_data_point("nIgG SBML", mcs * days_to_mcs, self.sbml.FullModel['nIgG'])
        self.plot_win7.add_data_point("sIgM SBML", mcs * days_to_mcs, self.sbml.FullModel['sIgM'])
        self.plot_win7.add_data_point("sIgG SBML", mcs * days_to_mcs, self.sbml.FullModel['sIgG'])

        if plot_Epithelial_cells:
            self.plot_win.add_data_point("ODEE", mcs * days_to_mcs,
                                         self.sbml.FullModel['E'] / self.sbml.FullModel['E0'])
            self.plot_win.add_data_point("ODEEv", mcs * days_to_mcs,
                                         self.sbml.FullModel['Ev'] / self.sbml.FullModel['E0'])
            self.plot_win.add_data_point("ODED", mcs * days_to_mcs,
                                         self.sbml.FullModel['D'] / self.sbml.FullModel['E0'])
            self.plot_win.add_data_point("CC3DE", mcs * days_to_mcs,
                                         len(self.cell_list_by_type(self.E)) / self.initial_uninfected)
            self.plot_win.add_data_point("CC3DEv", mcs * days_to_mcs,
                                         len(self.cell_list_by_type(self.EV)) / self.initial_uninfected)
            self.plot_win.add_data_point("CC3DD", mcs * days_to_mcs,
                                         len(self.cell_list_by_type(self.D)) / self.initial_uninfected)
        if plot_Virus:
            secretor = self.get_field_secretor("Virus")
            self.field_virus = 0.0
            for cell in self.cell_list:
                self.field_virus += secretor.amountSeenByCell(cell)
            self.plot_win2.add_data_point("ODEV", mcs * days_to_mcs,
                                          self.sbml.FullModel['V'])
            self.plot_win2.add_data_point("CC3DV", mcs * days_to_mcs,
                                          self.field_virus)

        if plot_Tcells:
            self.plot_win3.add_data_point("CC3DTc", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.TCELL)))
            self.plot_win3.add_data_point("ODETc", mcs * days_to_mcs,
                                          self.sbml.FullModel['g'] * self.sbml.FullModel['Tc']
                                          / self.sbml.FullModel['E0'] * self.initial_uninfected)

        if plot_APC:
            self.plot_win4.add_data_point("ODEAPC", mcs * days_to_mcs,
                                          self.sbml.FullModel['Da'] / self.sbml.FullModel['E0'])
            self.plot_win4.add_data_point("CC3DAPC0", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.APC)) / self.initial_uninfected)
            self.plot_win4.add_data_point("CC3DAPC", mcs * days_to_mcs,
                                          self.shared_steppable_vars['Active_APCs'] / self.initial_uninfected)

        if plot_Dm:
            self.plot_win5.add_data_point("ODEAPC", mcs * days_to_mcs,
                                          self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'])

            # self.plot_win5.add_data_point("ODEAPC", mcs * days_to_mcs,
            #                               self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'])
            self.plot_win5.add_data_point("CC3DAPC", mcs * days_to_mcs,
                                          self.sbml.LymphModel['Dm'] / self.sbml.FullModel['E0'])

        if contact_cytotoxicity and plot_Contact_Cytotoxicity:
            self.plot_win6.add_data_point("ODECC", mcs * days_to_mcs, self.sbml.FullModel['J6'])
            # self.plot_win6.add_data_point("ODECC", mcs * days_to_mcs, self.shared_steppable_vars['ODE_Killing'])
            self.plot_win6.add_data_point("CC3DCC", mcs * days_to_mcs, self.shared_steppable_vars['Contact_Killing'])
