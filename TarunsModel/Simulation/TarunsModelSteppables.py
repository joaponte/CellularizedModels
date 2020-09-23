from cc3d.core.PySteppables import *
import numpy as np

min_to_mcs = 60.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs
days_to_simulate = 30.0

virus_infection_feedback = 3
use_LymphModel_outputs = True

contact_cytotoxicity = True

plot_Epithelial_cells = True
plot_Virus = True
plot_Tcells = True
plot_APC = True
plot_Dm = True
plot_Contact_Cytotoxicity = True
plot_viral_antibodies = True
plot_infected_antibodies = True
plot_current_virus_decay = True

# TODO: add graphs for all state variables


# changed J6, added drt = 30;

full_model_string = """
/////////////////////// Reactions //////////////////////////////

//Tissue infection

J1: D -> E; dE*E0; // E0 => initial number of E-cells
J2: E -> D; dE*E;
J3: E -> Ev; bE*V*E;
J4: Ev -> E; aE*Ev;
J5: Ev -> D; dEv*Ev;
J6: Ev -> D; kE*Tct*Ev/(KEv + Ev + Tct);

J7: Tc -> Tct; g*Tc;
J8: Tct ->; dC*Tct;
J9: Tct ->; dT1*Tct*Da/(Da+dT2);

J10: -> V; pV*Ev;
J11: V ->; cV*V;
J12: -> Da; bD*V*(D0-Da); // D0 => initial number of inactive D-cells
J13: Da ->; dD*Da;

// Systemic reactions in Lymph

J14: -> Dm; kD*Da;
J15: Dm ->; dDm*Dm;
J16: -> Tc; dC*Tc0; // Tc0 => initial number of naive Tc-cells
J17: Tc ->; dC*Tc;
J18:  -> Tc; rT1*Dm*Tc/(Dm+rT2); 
J19: Tc ->; dT1*Tc*Dm/(Dm+dT2);

J20: -> Th1; sTh1*Th1/(1+Th2)^2;
J21:  -> Th1; pTh1*Dm*(Th1^2)/(1+Th2)^2;
J22: Th1 ->; dTh1*Dm*(Th1^3)/(1+Th2);
J23: Th1 ->; mTh*Th1;
J24: -> Th2; sTh2*Th2/(1+Th2);
J25:  -> Th2; pTh2*(ro+Th1)*Dm*(Th2^2)/((1+Th2)*(1+Th1+Th2)) 
J26: Th2 ->; mTh*Th2;

// Antibody production

J27: -> B; dB*B0; // B0 => initial number of inactive B-cells
J28: B ->; dB*B;
J29:  -> B; rB1*B*(Dm+h*Th2)/(Dm+h*Th2+rB2);
J30: B -> Pss; pS*B;
J31: B -> Psn; pS*B;
J32: B -> Pls; pL*B*Th2;
J33: B -> Pln; pL*B*Th2;
J34: Pss ->; dS*Pss;
J35: Psn ->; dS*Psn;
J36: Pls ->; dL*Pls;
J37: Pln ->; dL*Pln;
J38: Pss -> Pls; d*(1-v)*Pss;
J39: Psn -> Pln; d*(1-v)*Psn;
J40:  -> Pss; b*v*Pss;
J41:  -> Pls; b*v*Pls;
J42:  -> Psn; b*v*Psn;
J43:  -> Pln; b*v*Pln;
J44: Pls -> Pss; d*(1-v)*Pls;
J45: Pln -> Psn; d*(1-v)*Pln; 

J46:  -> sIgM; pAS*Pss;
J47:  -> nIgM; pAS*Psn;
J48:  -> sIgG; pAS*Pls;
J49:  -> nIgG; pAS*Pln;
J50: sIgM ->; dM*sIgM;
J51: sIgG ->; dG*sIgG;
J52: nIgM ->; dM*nIgM;
J53: nIgG ->; dG*nIgG;

// Antibody feedback to tissue

J54: Ev + nIgM -> D; eE*Ev*nIgM; 
J55: Ev + nIgG -> D; eE*Ev*nIgG;
J56: V + sIgM ->; eV*V*sIgM;
J57: V + sIgG ->; eV*V*sIgG;


/////////////////// Parameters ///////////////////////

// Epithelial infection

dE=10^-3;
E0 = 5.0*10^5;
bE=7.0*10^-6;
dEv=0.12;
aE=5.0*10^-1;
pV=1.9;
cV=1.0;

// Dendritic cell infection, activation, migration

D0=10^4;
bD=10^-6;
dD=2.9;
kD = 0.5;
dDm = 0.5;

// Cytotoxic T-cell activation, proliferation

dC=2.0*10^-3;
Tc0=2.0*10^3;
rT1=3.5;
rT2=2.0*10^3;
dT1=1.0;
dT2=1.0;

// T-cell mediated Cytotoxicity

kE=1.19*10^-2 / 900;
g=0.15 * 900; // 900 rescaling for non-0 T cell pop in tissue
KEv = 500.0;

// Helper T-cell activation, proliferation

sTh1=1.0;
pTh1=0.012;
dTh1=0.001;
KTh1 =500.0;
mTh=0.0225;
sTh2=0.04;
pTh2=0.003;
ro=1.0;

// B-cell activation, proliferation, differentiation

dB=0.02;
B0=2.0*10^1;
rB1=4.5;
rB2=1*10^4;
h=1.0;

// Plasma cell proliferation, differentiation and antibody production

pS=3.0*10^-1;
pL=1.5*10^-4;
dS=0.2;
dL=0.02;
b=2.4*10^-2;
d=2.4*10^-2;
pAS=0.8*10^2;
pAL=1.2*10^2;
dG=0.5;
dM=2.0;

// Switching functions of the Plasma Cells

u =  0.5;
v =  0.5;

// Antibody activity: virus and cell killing

eE=0.0001;
eV=0.00018;

////////////////// Initial Conditions /////////////////////

E = 5.0*10^5; // Uninfected epithelial cells
Ev = 10.0; // Virus-infected epithelial cells
V = 1000.0; // 
Da = 0.0; // Infected-activated dendtritic cells in Tissue

Dm = 0.0; // Migrated dendritic cells in Lymph
Tc = 1.0; // Effector cytotoxic T-cells in Lymph
Tct = 0.0; // Effector cytotoxic T-cells in Tissue
Th1 = 1.0; // Type I helper T-cells
Th2 = 1.0; // Type II helper T-cells

B = 1.0; // Activated B-cells
pSs = 0.0; // SP-RBD-specific Short-living plasma cells
pLs = 0.0; // SP-RBD-specific Long-living plasma cells
pSn = 0.0; // NP-specific Short-living plasma cells
pLn = 0.0; // NP-specific Long-living plasma cells

sIgM = 0.0; // SP-RBD-specific IgM
sIgG = 0.0; // SP-RBD-specific IgG
nIgM = 0.0; // NP-specific IgM
nIgG = 0.0; // NP-specific IgG
"""

lymph_node_string = '''
// Systemic reactions in Lymph

J7: Tc -> Tct; g*Tc;
J8: Tct ->; dC*Tct;
J9: Tct ->; dT1*Tct*Da/(Da+dT2);

J14: -> Dm; 0.0*kD*Da; // MUST STAY SHUT OFF, DM IS AN INPUT Dm are apc in lymph
J15: Dm ->; dDm*Dm;
J16: -> Tc; dC*Tc0; // Tc0 => initial number of naive Tc-cells
J17: Tc ->; dC*Tc;
J18:  -> Tc; rT1*Dm*Tc/(Dm+rT2); 
J19: Tc ->; dT1*Tc*Dm/(Dm+dT2);

J20: -> Th1; sTh1*Th1/(1+Th2)^2;
J21:  -> Th1; pTh1*Dm*(Th1^2)/(1+Th2)^2;
J22: Th1 ->; dTh1*Dm*(Th1^3)/(1+Th2);
J23: Th1 ->; mTh*Th1;
J24: -> Th2; sTh2*Th2/(1+Th2);
J25:  -> Th2; pTh2*(ro+Th1)*Dm*(Th2^2)/((1+Th2)*(1+Th1+Th2)) 
J26: Th2 ->; mTh*Th2;

// Antibody production

J27: -> B; dB*B0; // B0 => initial number of inactive B-cells
J28: B ->; dB*B;
J29:  -> B; rB1*B*(Dm+h*Th2)/(Dm+h*Th2+rB2);
J30: B -> Pss; pS*B;
J31: B -> Psn; pS*B;
J32: B -> Pls; pL*B*Th2;
J33: B -> Pln; pL*B*Th2;
J34: Pss ->; dS*Pss;
J35: Psn ->; dS*Psn;
J36: Pls ->; dL*Pls;
J37: Pln ->; dL*Pln;
J38: Pss -> Pls; d*(1-v)*Pss;
J39: Psn -> Pln; d*(1-v)*Psn;
J40:  -> Pss; b*v*Pss;
J41:  -> Pls; b*v*Pls;
J42:  -> Psn; b*v*Psn;
J43:  -> Pln; b*v*Pln;
J44: Pls -> Pss; d*(1-v)*Pls;
J45: Pln -> Psn; d*(1-v)*Pln; 

J46:  -> sIgM; pAS*Pss;
J47:  -> nIgM; pAS*Psn;
J48:  -> sIgG; pAS*Pls;
J49:  -> nIgG; pAS*Pln;
J50: sIgM ->; dM*sIgM;
J51: sIgG ->; dG*sIgG;
J52: nIgM ->; dM*nIgM;
J53: nIgG ->; dG*nIgG;

// Antibody feedback to tissue

J54: Ev + nIgM -> D; eE*Ev*nIgM; 
J55: Ev + nIgG -> D; eE*Ev*nIgG;
J56: V + sIgM ->; eV*V*sIgM;
J57: V + sIgG ->; eV*V*sIgG;


/////////////////// Parameters ///////////////////////

// Epithelial infection

dE=10^-3;
E0 = 5.0*10^5;
bE=7.0*10^-6;
dEv=0.12;
aE=5.0*10^-1;
pV=1.9;
cV=1.0;

// Dendritic cell infection, activation, migration

D0=10^4;
bD=10^-6;
dD=2.9;
kD = 0.5;
dDm = 0.5;

// Cytotoxic T-cell activation, proliferation

dC=2.0*10^-3;
Tc0=2.0*10^3;
rT1=3.5;
rT2=2.0*10^3;
dT1=1.0;
dT2=1.0;

// T-cell mediated Cytotoxicity

kE=1.19*10^-2 / 900;
g=0.15 * 900; // 900 rescaling for non-0 T cell pop in tissue
KEv = 500.0;

// Helper T-cell activation, proliferation

sTh1=1.0;
pTh1=0.012;
dTh1=0.001;
KTh1 =500.0;
mTh=0.0225;
sTh2=0.04;
pTh2=0.003;
ro=1.0;

// B-cell activation, proliferation, differentiation

dB=0.02;
B0=2.0*10^1;
rB1=4.5;
rB2=1*10^4;
h=1.0;

// Plasma cell proliferation, differentiation and antibody production

pS=3.0*10^-1;
pL=1.5*10^-4;
dS=0.2;
dL=0.02;
b=2.4*10^-2;
d=2.4*10^-2;
pAS=0.8*10^2;
pAL=1.2*10^2;
dG=0.5;
dM=2.0;

// Switching functions of the Plasma Cells

u =  0.5;
v =  0.5;

// Antibody activity: virus and cell killing

eE=0.0001;
eV=0.00018;

////////////////// Initial Conditions /////////////////////

E = 5.0*10^5; // Uninfected epithelial cells
Ev = 10.0; // Virus-infected epithelial cells
V = 1000.0; // 
Da = 0.0; // Infected-activated dendtritic cells in Tissue

Dm = 0.0; // Migrated dendritic cells in Lymph
Tc = 1.0; // Effector cytotoxic T-cells in Lymph
Tct = 0.0; // Effector cytotoxic T-cells in Tissue
Th1 = 1.0; // Type I helper T-cells
Th2 = 1.0; // Type II helper T-cells

B = 1.0; // Activated B-cells
pSs = 0.0; // SP-RBD-specific Short-living plasma cells
pLs = 0.0; // SP-RBD-specific Long-living plasma cells
pSn = 0.0; // NP-specific Short-living plasma cells
pLn = 0.0; // NP-specific Long-living plasma cells

sIgM = 0.0; // SP-RBD-specific IgM
sIgG = 0.0; // SP-RBD-specific IgG
nIgM = 0.0; // NP-specific IgM
nIgG = 0.0; // NP-specific IgG
'''


class TarunsModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def get_ode_parameters(self):

        ode_params = {'dE': 1e-3, 'E0': 5e-2, 'bE': 3e-6, 'aE': 5e-2, 'V0': 10.0, 'pV': 19.0, 'cV': 1.0,
                      'kE': 1.19e-2 / 900, 'g': 0.15 * 900,  # rescaling for non - 0 T cell pop in tissue (/900*900)
                      'tC': 0.5,  # not included in the model
                      'eE': 0.05, 'eV': 16.0, 'D0': 1e3, 'bD': 1e-7,
                      'dD': 2.9, 'kD': 0.5,
                      'tD': 1.0,  # not included in the model
                      'dDm': 0.5, 'dC': 10.1e-3, 'Tc0': 5e2, 'rT1': 1.3, 'rT2': 100.0,
                      'dT1': 5.0, 'dT2': 200.0, 'sTh1': 0.25, 'pTh1': 0.4, 'dTh1': 0.03, 'mTh': 0.25, 'sTh2': 0.001,
                      'pTh2': 0.0022, 'ro': 0.2, 'dB': 0.0009, 'B0': 1e3, 'rB1': 1e2, 'rB2': 2e5, 'h': 100, 'pS': 0.1,
                      'pL': 8e-3, 'dS': 0.002, 'dL': 0.02, 'b': 2.4e-4, 'd': 2.4e-2, 'pAS': 0.8, 'pAL': 1.2, 'dG': 0.04,
                      'dM': 0.2, 'pT2': 600.0, 'v': 0.5,
                      'proport_init_inf': 0.01  # unused
                      }
        ode_params['E'] = ode_params['E0'] - 1
        ode_params['Ev'] = 1
        ode_params['V'] = 0
        ode_params['Th1'] = 10.0
        ode_params['Th2'] = 10.0
        ode_params['drt'] = 3.0

        return ode_params

    def decay_exp_prob(self, p):
        return 1 - np.exp(-p)

    def seed_epithelial_sheet(self):
        for x in range(0, self.dim.x, int(3)):
            for y in range(0, self.dim.y, int(3)):
                cell = self.new_cell(self.E)
                self.cellField[x:x + int(3), y:y + int(3), 0] = cell

    def init_variables(self):
        # Updating max simulation steps using scaling factor to simulate 10 days
        self.add_free_floating_antimony(model_string=full_model_string, model_name='FullModel', step_size=days_to_mcs)
        self.add_free_floating_antimony(model_string=lymph_node_string, model_name='LymphModel', step_size=days_to_mcs)
        self.get_xml_element('simulation_steps').cdata = days_to_simulate / days_to_mcs
        # self.initial_uninfected = len(self.cell_list_by_type(self.E))
        self.get_xml_element('virus_decay').cdata = self.sbml.FullModel['cV'] * days_to_mcs
        lattice_factor = .14
        self.max_virus_gamma = .75 * float(self.get_xml_element("virus_D").cdata) / lattice_factor
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
                cell.targetVolume = cell.volume + .5
                cell.lambdaVolume = cell.volume * 3
                cell.dict['Activation_State'] = 0  # in tissue
                numberAPC += 1

    def new_T_cell_in_location(self, x, y, z):
        # print('in net t cell func')
        cell = self.new_cell(self.TCELL)
        cell.dict['body_count'] = 0
        self.cell_field[x:x + 3, y:y + 3, z] = cell
        cell.targetVolume = cell.volume + .5
        cell.lambdaVolume = cell.volume * 3
        cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
        cd.setLambda(0)
        cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
        return cell

    def start(self):
        self.seed_epithelial_sheet()

        self.init_variables()
        p_Ev = self.sbml.FullModel['Ev'] / self.sbml.FullModel['Ev']
        # print('@@@@@@@@@@@@@@@@@@@@@@@\n', p_Ev, '\n@@@@@@@@@@@@@@@@@@@@@@')
        proport_inf = False
        if proport_inf:
            for cell in self.cell_list_by_type(self.E):
                if np.random.random() <= p_Ev:
                    # print('@@@@@@@@@@@@@@@@@@@@@@@\ninit infected\n@@@@@@@@@@@@@@@@@@@@@@')
                    cell.type = self.EV
        else:
            cells = [cell for cell in self.cell_list_by_type(self.E)]
            np.random.shuffle(cells)
            cell = cells[0]
            cell.type = self.EV
        self.initial_uninfected = len(self.cell_list_by_type(self.E))
        self.seed_APC_cells()

    def J1_DtoE(self, cell):
        ## Transition from D to E
        # J1: D -> E; dE*E0; // E0 => initial number of E-cells

        # 1 - np.exp(-self.sbml.FullModel['dE'] * days_to_mcs) ~ 1 - ( 1 - self.sbml.FullModel['dE'] * days_to_mcs)
        # = self.sbml.FullModel['dE'] * days_to_mcs
        p = self.decay_exp_prob(self.sbml.FullModel['dE'] * self.sbml.FullModel['E0'] * days_to_mcs)
        if p > np.random.random():
            cell.type = self.E

    def J2_EtoD(self, cell):
        ## Transition from E to D
        # J2: E -> D; dE * E;
        # print('J2',self.sbml.FullModel['dE'] * days_to_mcs)
        if self.decay_exp_prob(self.sbml.FullModel['dE'] * days_to_mcs) > np.random.random():
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
        p_EtoEv = self.decay_exp_prob(bE * V)
        # p_EtoEv = bE * V

        if p_EtoEv > np.random.random():
            cell.type = self.EV
            # print('infected', cell.id)
            # print('prob infect', p_EtoEv)
            # print('old prob inf', bE * V)

    def J4_EvtoE(self, cell):
        ## Transition from Ev to E
        # J4: Ev -> E; aE*Ev;
        if self.decay_exp_prob(self.sbml.FullModel['aE'] * days_to_mcs) > np.random.random():
            cell.type = self.E
            # print('became healthy', cell.id)

    def J5_EvtoD(self, cell):
        ## Transition from Ev to D
        # J5: Ev -> D; dEv*Ev;
        dE = self.sbml.FullModel['dEv'] * days_to_mcs
        p_EvtoD = self.decay_exp_prob(dE)
        if p_EvtoD > np.random.random():
            cell.type = self.D
            # print('j5', p_EvtoD)
            # print('died from j5', cell.id)

    def J6_EvtoD_ODE_killing(self, cell, Tct):
        ## Transition from Ev to D
        # Tct is new tc in tissue, was g*Tc
        # J6: Ev -> D; kE*Tct*Ev/(KEv + Ev + Tct);
        p = self.sbml.FullModel['kE'] * Tct * 1 / (self.sbml.FullModel['KEv'] + 1 + Tct)
        p_EvtoD = self.decay_exp_prob(p)
        if p_EvtoD > np.random.random():
            if not contact_cytotoxicity:
                cell.type = self.D
            self.shared_steppable_vars['ODE_Killing'] += 1

    def J7_Tct_seeding(self):
        # J7: -> Tct; g*Tc;
        use_LymphModel_outputs = False
        if use_LymphModel_outputs:
            n_tct = self.sbml.LymphModel['Tct'] * \
                    self.initial_uninfected / self.sbml.FullModel['E0']
        else:
            n_tct = self.sbml.FullModel['Tct'] * \
                    self.initial_uninfected / self.sbml.FullModel['E0']
        lm_n_tct = self.sbml.LymphModel['Tct'] * \
                    self.initial_uninfected / self.sbml.FullModel['E0']
        print(n_tct, len(self.cell_list_by_type(self.TCELL)), lm_n_tct)
        if int(n_tct) > len(self.cell_list_by_type(self.TCELL)):
            for i in range(int(round(n_tct - len(self.cell_list_by_type(self.TCELL))))):
                cell = False
                while not cell:
                    x = np.random.randint(0, self.dim.x - 3)
                    y = np.random.randint(0, self.dim.y - 3)
                    if not self.cell_field[x, y, 1]:
                        cell = self.new_T_cell_in_location(x, y, 1)
        return

    def J8_J9_Tct_death(self, cell, tissue_apc):
        # change if to be on the call to this function

        # J8: Tct ->; dC * Tct;
        # J9: Tct ->; dT1 * Tct * Da / (Da + dT2);

        dC = self.sbml.FullModel['dC'] * days_to_mcs
        dT1 = self.sbml.FullModel['dT1'] * days_to_mcs
        dT2 = self.sbml.FullModel['dT2'] * days_to_mcs * self.initial_uninfected / self.sbml.FullModel['E0']
        p = dC * 1 + dT1 * 1 * tissue_apc / (dT2 + tissue_apc)
        p = self.decay_exp_prob(p)
        if p > np.random.random():
            cell.targetVolume = 0
            cell.lambdaVolume = 99
        return

    def J10_virus_production(self):
        ## Virus Production
        # J7: -> V; pV*Ev;
        virus_production = 0.0
        pV = self.sbml.FullModel['pV'] * days_to_mcs * self.sbml.FullModel['E0'] / self.initial_uninfected
        for cell in self.cell_list_by_type(self.EV):
            release = self.virus_secretor.secreteInsideCellTotalCount(cell, pV / cell.volume)
            virus_production += abs(release.tot_amount)
        return virus_production

    def old_J8_virus_decay(self):
        # no longer used, now update_scalar_virus does all the work
        ## Virus Decay
        # J8: V ->; cV*V;
        cV = self.sbml.FullModel['cV'] * days_to_mcs
        virus_decay = cV * self.scalar_virus

    def J12_J13_APC_activation_deactivation(self):
        tissue_apc = 0
        lymph_apc = 0
        node_apc = 0
        just_moved_in_to_node = 0
        # activated_APC_count = 0
        for cell in self.cell_list_by_type(self.APC):
            # if cell.dict == 0 -> in tissue; keep logic (J12: -> Da; bD*V*(D0-Da);) but replace cell.dict[
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
                # J12: -> Da; bD*V*(D0-Da);
                tissue_apc += 1
                bD = (self.sbml.FullModel['bD'] * days_to_mcs) * (self.sbml.FullModel['D0'] * self.initial_uninfected /
                                                                  self.sbml.FullModel['E0'])
                # V should be local instead of the total virus
                V = self.virus_secretor.amountSeenByCell(cell) * self.initial_uninfected
                # V = self.sbml.FullModel['V'] / self.sbml.FullModel['E0'] * self.initial_uninfected

                p_tissue_2_lymph = self.decay_exp_prob(bD * V)
                if p_tissue_2_lymph > np.random.random():
                    cell.dict['Activation_State'] = 1
                    tissue_apc -= 1
                    lymph_apc += 1
            elif cell.dict['Activation_State'] == 1:  # in lymph
                # transport of apc to lymph
                lymph_apc += 1
                # todo: rethink this probability and find a good parameter for transport
                p_lymph_2_node = 1 / 50  # PLACEHOLDER
                if p_lymph_2_node > np.random.random():
                    cell.dict['Activation_State'] = 2
                    lymph_apc -= 1
                    node_apc += 1
                    just_moved_in_to_node += 1
            elif cell.dict['Activation_State'] == 2:
                # J13: Da ->; dD*Da;
                node_apc += 1
                dD = self.sbml.FullModel['dD'] * days_to_mcs
                p_lymph_2_tissue = self.decay_exp_prob(dD)
                if p_lymph_2_tissue > np.random.random():
                    cell.dict['Activation_State'] = 0
                    node_apc -= 1
                    tissue_apc += 1
        self.shared_steppable_vars['Active_APCs'] = node_apc  # lymph_apc + node_apc
        return tissue_apc, lymph_apc, node_apc, just_moved_in_to_node

    def J14_APC_travel_Lymph_Model_input(self, just_moved_in_to_node):
        ## APC "travel" to lymph node
        # we'll implement it as a signal that is proportional to the activated apcs
        # J14: -> Dm; kD*Da; // Dm are apc in lymph
        # print(node_apc)
        self.sbml.LymphModel['Dm'] += max(0, (just_moved_in_to_node * self.sbml.FullModel['D0'] *
                                              self.initial_uninfected / self.sbml.FullModel['E0']))

    def old_J13_Tcell_stable_population_seeding(self):
        ## Tcell seeding
        # J13: -> Tc; dc*g*Tc0;
        if use_LymphModel_outputs:
            dc = self.sbml.LymphModel['dC'] * days_to_mcs
            g = self.sbml.LymphModel['g']
            Tc0 = self.sbml.LymphModel['Tc0'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        else:
            dc = self.sbml.FullModel['dC'] * days_to_mcs
            g = self.sbml.FullModel['g']
            Tc0 = self.sbml.FullModel['Tc0'] / self.sbml.FullModel['E0'] * self.initial_uninfected
        cell = False
        if self.decay_exp_prob(dc * g * Tc0) > np.random.random():
            while not cell:
                x = np.random.randint(0, self.dim.x - 3)
                y = np.random.randint(0, self.dim.y - 3)
                if not self.cell_field[x, y, 1]:
                    cell = self.new_T_cell_in_location(x, y, 1)

    def old_J14_Tcell_clearance(self):
        ## Clearance of Tcells
        # J14: Tc ->; dc*Tc;

        for cell in self.cell_list_by_type(self.TCELL):
            dc = self.decay_exp_prob(self.sbml.FullModel['dC'] * days_to_mcs)
            if dc > np.random.random():
                cell.targetVolume = 0.0

    def old_J15a_Tcell_inflamatory_seeding(self):  # REMOVE!!!
        # equation changed no longer used
        return
        ## Tcell seeding
        # J15a: Dm -> Tc; Dm ;
        # J15a:  -> Tc; Dm ;
        g = self.sbml.FullModel['g']

        if use_LymphModel_outputs:
            Dm = self.sbml.LymphModel['Dm'] * days_to_mcs / self.sbml.FullModel['E0'] * self.initial_uninfected
        else:
            Dm = self.sbml.FullModel['Dm'] * days_to_mcs / self.sbml.FullModel['E0'] * self.initial_uninfected

        if self.decay_exp_prob(Dm * g) > np.random.random():
            cells_to_seed = max(1, round(Dm * g))
            for i in range(cells_to_seed):
                cell = False
                while not cell:
                    x = np.random.randint(0, self.dim.x - 3)
                    y = np.random.randint(0, self.dim.y - 3)
                    if not self.cell_field[x, y, 1]:
                        cell = self.new_T_cell_in_location(x, y, 1)

    def old_J15b_Tcell_inflamatory_seeding(self):
        ## Tcell seeding
        # J15b:  -> Tc; rT1 * Dm * Tc / (Dm + pT2)
        if use_LymphModel_outputs:
            rT1 = self.sbml.LymphModel['rT1'] * days_to_mcs
            Dm = self.sbml.LymphModel['Dm'] / self.sbml.FullModel['E0'] * self.initial_uninfected
            Tc = self.sbml.LymphModel['Tc'] / self.sbml.FullModel['E0'] * self.initial_uninfected
            pT2 = self.sbml.LymphModel['pT2'] / self.sbml.FullModel['E0'] * self.initial_uninfected
            g = self.sbml.LymphModel['g']
        else:
            rT1 = self.sbml.FullModel['rT1'] * days_to_mcs
            Dm = self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'] * self.initial_uninfected
            Tc = self.sbml.FullModel['Tc'] / self.sbml.FullModel['E0'] * self.initial_uninfected
            pT2 = self.sbml.FullModel['pT2'] / self.sbml.FullModel['E0'] * self.initial_uninfected
            g = self.sbml.FullModel['g']

        if self.decay_exp_prob(g * rT1 * Dm * Tc / (Dm + pT2)) > np.random.random():
            cells_to_seed = max(1, round(g * rT1 * Dm * Tc / (Dm + pT2)))
            for i in range(cells_to_seed):
                cell = False
                while not cell:
                    x = np.random.randint(0, self.dim.x - 3)
                    y = np.random.randint(0, self.dim.y - 3)
                    if not self.cell_field[x, y, 1]:
                        cell = self.new_T_cell_in_location(x, y, 1)

    def old_J15_Tcell_inflamatory_seeding(self):
        # self.J15a_Tcell_inflamatory_seeding()
        self.old_J15b_Tcell_inflamatory_seeding()

    def old_J16_Tcell_clearance(self):
        # TODO: NEEDs RESCALING
        # Tcell clearance
        # J16: Tc ->; dT1 * Tc * Ev/(Ev+dT2)
        Ev = len(self.cell_list_by_type(self.EV))
        dT1 = self.sbml.FullModel['dT1'] * days_to_mcs * 3
        dT2 = self.sbml.FullModel['dT1'] * days_to_mcs / self.sbml.FullModel['E0'] * self.initial_uninfected
        for cell in self.cell_list_by_type(self.TCELL):
            if self.decay_exp_prob(dT1 * Ev / (Ev + dT2)) > np.random.random():
                cell.targetVolume = 0.0

    def old_J49_Ev2D_from_nIgM(self, cell):
        ## Transition from Ev to D
        # J49: Ev -> D;   eE * Ev * nIgM;
        if use_LymphModel_outputs:
            nIgM = self.sbml.LymphModel['nIgM']
            eE = self.sbml.LymphModel['eE']
        else:
            nIgM = self.sbml.FullModel['nIgM']
            eE = self.sbml.FullModel['eE']
        p_Ev2D = nIgM * eE * days_to_mcs  # * self.initial_uninfected / self.sbml.FullModel['E0']
        p_Ev2D = (p_Ev2D)
        if p_Ev2D > np.random.random():
            cell.type = self.D
            # print('killed by anti bodies M', cell.id)
            # print(p_Ev2D)
        return

    def old_J50_Ev2D_from_nIgG(self, cell):
        # J50: Ev -> D;    eE * Ev * nIgG;
        if use_LymphModel_outputs:
            nIgG = self.sbml.LymphModel['nIgG']
            eE = self.sbml.LymphModel['eE']
        else:
            nIgG = self.sbml.FullModel['nIgG']
            eE = self.sbml.FullModel['eE']
        p_Ev2D = nIgG * eE * days_to_mcs  # * self.initial_uninfected / self.sbml.FullModel['E0']
        p_Ev2D = self.decay_exp_prob(p_Ev2D)
        if p_Ev2D > np.random.random():
            cell.type = self.D
            # print('killed by anti bodies G', cell.id)
            # print(p_Ev2D)
        return

    def J54_J55_antibody_cell_death(self):
        # J54: Ev + nIgM -> D;    eE * Ev * nIgM;
        # J55: Ev + nIgG -> D;    eE * Ev * nIgG;
        if use_LymphModel_outputs:
            nIgM = self.sbml.LymphModel['nIgM']
            nIgG = self.sbml.LymphModel['nIgG']
            eE = self.sbml.LymphModel['eE']
        else:
            nIgM = self.sbml.FullModel['nIgM']
            nIgG = self.sbml.FullModel['nIgG']
            eE = self.sbml.FullModel['eE']
        p = days_to_mcs * eE * (nIgG + nIgM)
        p = self.decay_exp_prob(p)
        for cell in self.cell_list_by_type(self.EV):
            if p > np.random.random():
                cell.type = self.D
        return

    def J56_J57_update_virus_decay(self):
        # J8: V ->; cV*V;
        # J56: V + sIgM ->; eV*V*sIgM;
        # J57: V + sIgG ->; eV*V*sIgG;

        if use_LymphModel_outputs:
            gamma = (self.sbml.LymphModel['eV'] * (self.sbml.LymphModel['sIgM'] + self.sbml.LymphModel['nIgM']) +
                     self.sbml.LymphModel['cV']) * days_to_mcs
        else:
            gamma = (self.sbml.FullModel['eV'] * (self.sbml.FullModel['sIgM'] + self.sbml.FullModel['nIgM']) +
                     self.sbml.FullModel['cV']) * days_to_mcs

        effective_gamma = min(self.max_virus_gamma, gamma)
        # effective_gamma = (gamma_sIgM + gamma_sIgG + self.sbml.FullModel['cV']) * days_to_mcs
        self.get_xml_element('virus_decay').cdata = effective_gamma

        return effective_gamma

    def update_scalar_virus(self, virus_production, current_gamma):
        # does the work that J8 used to do
        if self.mcs % 50 == 0 and self.mcs > 0:
            # do the actual integral
            self.scalar_virus = self.virus_secretor.totalFieldIntegral()
            return
        virus_decay = current_gamma * self.scalar_virus
        self.scalar_virus += virus_production - virus_decay
        self.shared_steppable_vars['scalar_virus'] = self.scalar_virus
        return

    def lymph_model_input(self, Ev, Da):
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
                        cell.dict['body_count'] += 1
                        # print('killed by contact killing', cell.id)

    def test_prints(self):
        print('taruns\tcc3d')
        print('Tc', self.sbml.FullModel['Tc'], self.sbml.LymphModel['Tc'])
        print('Th1', self.sbml.FullModel['Th1'], self.sbml.LymphModel['Th1'])
        print('Th2', self.sbml.FullModel['Th2'], self.sbml.LymphModel['Th2'])
        # print('Da', self.sbml.FullModel['Da'], self.sbml.LymphModel['Da'])
        print('Dm', self.sbml.FullModel['Dm'], self.sbml.LymphModel['Dm'])

    def step(self, mcs):
        # self.test_prints()
        # print(mcs)
        self.shared_steppable_vars['Contact_Killing'] = 0

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

        self.J7_Tct_seeding()
        Tct = len(self.cell_list_by_type(self.TCELL)) / self.initial_uninfected * self.sbml.FullModel['E0']

        for cell in self.cell_list_by_type(self.EV):
            self.J6_EvtoD_ODE_killing(cell, Tct)

        # for cell in self.cell_list_by_type(self.EV):
        #     self.old_J49_Ev2D_from_nIgM(cell)
        #
        # for cell in self.cell_list_by_type(self.EV):
        #     self.old_J50_Ev2D_from_nIgG(cell)

        self.J54_J55_antibody_cell_death()

        virus_production = self.J10_virus_production()

        current_gamma = self.J56_J57_update_virus_decay()
        self.update_scalar_virus(virus_production, current_gamma)
        # activated_APC_count = 0

        tissue_apc, lymph_apc, node_apc, just_moved_in_to_node = self.J12_J13_APC_activation_deactivation()
        if use_LymphModel_outputs:
            Da = tissue_apc
        else:
            Da = self.sbml.FullModel['Da']
        for cell in self.cell_list_by_type(self.TCELL):
            self.J8_J9_Tct_death(cell, Da)

        self.J14_APC_travel_Lymph_Model_input(just_moved_in_to_node)

        # self.old_J13_Tcell_stable_population_seeding()
        #
        # self.old_J14_Tcell_clearance()
        #
        # self.old_J15_Tcell_inflamatory_seeding()
        #
        # self.old_J16_Tcell_clearance()

        ## Tcell Contact Killing
        for cell in self.cell_list_by_type(self.TCELL):
            self.contact_killing(cell)

        # update lymph model inputs
        self.lymph_model_input(len(self.cell_list_by_type(self.EV)), node_apc)

        # if use_LymphModel_outputs:
        #     self.lymph_model_input(len(self.cell_list_by_type(self.EV)), node_apc)
        # else:
        #     self.lymph_model_input(len(self.cell_list_by_type(self.EV)), self.sbml.FullModel['Da'])
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
        lambda_chemotaxis = 250.0 * 2
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
                                                     grid=True)

            self.plot_win.add_plot("ODEE", style='Lines', color='blue', size=5)
            self.plot_win.add_plot("ODEEv", style='Lines', color='yellow', size=5)
            self.plot_win.add_plot("ODED", style='Lines', color='red', size=5)

            self.plot_win.add_plot("CC3DE", style='Dots', color='blue', size=5)
            self.plot_win.add_plot("CC3DEv", style='Dots', color='yellow', size=5)
            self.plot_win.add_plot("CC3DD", style='Dots', color='red', size=5)

        if plot_Virus:
            self.plot_win2 = self.add_new_plot_window(title='Virus',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Virus', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True)

            self.plot_win2.add_plot("ODEV", style='Lines', color='yellow', size=5)
            self.plot_win2.add_plot("CC3DV", style='Dots', color='red', size=5)

        if plot_Tcells:
            self.plot_win3 = self.add_new_plot_window(title='Epithelial Tcells',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True)

            self.plot_win3.add_plot("ODETc", style='Lines', color='yellow', size=5)
            self.plot_win3.add_plot("CC3DTc", style='Dots', color='red', size=5)

        if plot_APC:
            self.plot_win4 = self.add_new_plot_window(title='Epithelial APC',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True)

            self.plot_win4.add_plot("ODEAPC", style='Lines', color='yellow', size=5)
            self.plot_win4.add_plot("CC3DAPC0", style='Dots', color='orange', size=5)
            self.plot_win4.add_plot("CC3DAPC", style='Dots', color='red', size=5)

        if plot_Dm:
            self.plot_win5 = self.add_new_plot_window(title='Lymph node APC',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True)

            self.plot_win5.add_plot("ODEAPC", style='Lines', color='yellow', size=5)
            # self.plot_win5.add_plot("CC3DAPC0", style='Lines', color='orange', size=5)

            self.plot_win5.add_plot("CC3DAPC", style='Dots', color='red', size=5)

        if contact_cytotoxicity and plot_Contact_Cytotoxicity:
            self.plot_win6 = self.add_new_plot_window(title='Contact Cytotoxicity',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Number of Cells', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True)

            self.plot_win6.add_plot("ODECC", style='Lines', color='yellow', size=5)
            self.plot_win6.add_plot("CC3DCC", style='Dots', color='red', size=5)
            self.plot_win6.add_plot("body_count", style='Dots', color='blue', size=5)
        if plot_viral_antibodies:
            self.plot_win7 = self.add_new_plot_window(title='Viral Anti-Bodies',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Conc', x_scale_type='linear',
                                                      y_scale_type='log',
                                                      grid=True)

            self.plot_win7.add_plot("sIgM control SBML", style='Lines', color='yellow', size=5)
            self.plot_win7.add_plot("sIgG control SBML", style='Lines', color='magenta', size=5)

            self.plot_win7.add_plot("sIgM CC3D SBML", style='Dots', color='yellow', size=5)
            self.plot_win7.add_plot("sIgG CC3D SBML", style='Dots', color='magenta', size=5)

        if plot_infected_antibodies:
            self.plot_win8 = self.add_new_plot_window(title='Infected Cell Anti-Bodies',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='Conc', x_scale_type='linear',
                                                      y_scale_type='log',
                                                      grid=True)
            self.plot_win8.add_plot("nIgM control SBML", style='Lines', color='yellow', size=5)
            self.plot_win8.add_plot("nIgG control SBML", style='Lines', color='magenta', size=5)

            self.plot_win8.add_plot("nIgM CC3D SBML", style='Dots', color='yellow', size=5)
            self.plot_win8.add_plot("nIgG CC3D SBML", style='Dots', color='magenta', size=5)
        if plot_current_virus_decay:
            self.plot_win9 = self.add_new_plot_window(title='effective virus decay',
                                                      x_axis_title='Time (days)',
                                                      y_axis_title='gamma', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True)
            self.plot_win9.add_plot('ODE effective decay', style='Lines', color='yellow', size=5)
            self.plot_win9.add_plot('CC3D effective decay', style='Dots', color='red', size=5)

    def step(self, mcs):

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
            for cell in self.cell_list_by_type(self.E, self.EV, self.D):
                self.field_virus += secretor.amountSeenByCell(cell)
            self.plot_win2.add_data_point("ODEV", mcs * days_to_mcs,
                                          self.sbml.FullModel['V'])
            # self.plot_win2.add_data_point("CC3DV", mcs * days_to_mcs, 2 * self.field_virus)
            self.plot_win2.add_data_point("CC3DV", mcs * days_to_mcs, self.shared_steppable_vars['scalar_virus'])

        if plot_Tcells:
            self.plot_win3.add_data_point("CC3DTc", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.TCELL)))
            self.plot_win3.add_data_point("ODETc", mcs * days_to_mcs, self.sbml.FullModel['Tct'] *
                                          self.initial_uninfected / self.sbml.FullModel['E0'])

        if plot_APC:
            self.plot_win4.add_data_point("ODEAPC", mcs * days_to_mcs,
                                          self.sbml.FullModel['Da'] / self.sbml.FullModel['E0'])
            self.plot_win4.add_data_point("CC3DAPC0", mcs * days_to_mcs,
                                          len(self.cell_list_by_type(self.APC)) / self.initial_uninfected)
            self.plot_win4.add_data_point("CC3DAPC", mcs * days_to_mcs,
                                          self.shared_steppable_vars['Active_APCs'] / self.initial_uninfected)
            # self.plot_win4.add_data_point("CC3DAPC", mcs * days_to_mcs,
            #                               self.sbml.LymphModel['Da'] / self.sbml.FullModel['E0'])

        if plot_Dm:
            self.plot_win5.add_data_point("ODEAPC", mcs * days_to_mcs,
                                          self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'])

            # self.plot_win5.add_data_point("ODEAPC", mcs * days_to_mcs,
            #                               self.sbml.FullModel['Dm'] / self.sbml.FullModel['E0'])
            self.plot_win5.add_data_point("CC3DAPC", mcs * days_to_mcs,
                                          self.sbml.LymphModel['Dm'] / self.sbml.FullModel['E0'])

        if contact_cytotoxicity and plot_Contact_Cytotoxicity:
            self.plot_win6.add_data_point("ODECC", mcs * days_to_mcs,
                                          self.sbml.FullModel['J6'] / self.sbml.FullModel['E0'])
            # g * Tc
            # self.plot_win6.add_data_point("ODECC", mcs * days_to_mcs, self.shared_steppable_vars['ODE_Killing'])
            self.plot_win6.add_data_point("CC3DCC", mcs * days_to_mcs,
                                          self.shared_steppable_vars['Contact_Killing'] / self.initial_uninfected)

        if plot_viral_antibodies:
            self.plot_win7.add_data_point("sIgM control SBML", mcs * days_to_mcs, self.sbml.FullModel['sIgM'])
            self.plot_win7.add_data_point("sIgG control SBML", mcs * days_to_mcs, self.sbml.FullModel['sIgG'])

            self.plot_win7.add_data_point("sIgM CC3D SBML", mcs * days_to_mcs, self.sbml.LymphModel['sIgM'])
            self.plot_win7.add_data_point("sIgG CC3D SBML", mcs * days_to_mcs, self.sbml.LymphModel['sIgG'])
        if plot_infected_antibodies:
            self.plot_win8.add_data_point("nIgM control SBML", mcs * days_to_mcs, self.sbml.FullModel['nIgM'])
            self.plot_win8.add_data_point("nIgG control SBML", mcs * days_to_mcs, self.sbml.FullModel['nIgG'])

            self.plot_win8.add_data_point("nIgM CC3D SBML", mcs * days_to_mcs, self.sbml.LymphModel['nIgM'])
            self.plot_win8.add_data_point("nIgG CC3D SBML", mcs * days_to_mcs, self.sbml.LymphModel['nIgG'])
        if plot_current_virus_decay:
            ode_g = (self.sbml.FullModel['eV'] * (self.sbml.FullModel['sIgM'] + self.sbml.FullModel['nIgM']) +
                     self.sbml.FullModel['cV']) * days_to_mcs
            if ode_g != 0:
                self.plot_win9.add_data_point('ODE effective decay', mcs * days_to_mcs, ode_g)
            cc3d_g = self.get_xml_element('virus_decay').cdata
            if cc3d_g != 0:
                self.plot_win9.add_data_point('CC3D effective decay', mcs * days_to_mcs, cc3d_g)
