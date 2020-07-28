import math

# Name of Antimony/SBML model of immune cell recruitment
ir_model_name = 'TcellModel'

# Name of Antimony/SBML model of immune cell recruitment in ODE form
ir_model_name_ode = 'TcellModelODE'

# Name of Antimony/SBML model of immune cell recruitment in original form
ir_model_name_orig = 'TcellModelOriginal'

# Key to reference of ImmuneRecruitmentSteppable instance in shared global dictionary
ir_steppable_key = 'ir_steppable'

# Key to reference of SimDataSteppable instance in shared global dictionary
simdata_steppable_key = 'simdata_steppable'


def immune_recruitment_model_string(num_ec=1E7, sites_per_cell=1, time_conv=1):
    model_string = f'''model {ir_model_name}()
    // Scaling coefficients
    s_v = {num_ec} / 1E7;
    s_l = 1 / 1E7 / {sites_per_cell};
    s_t = {time_conv};

    // Reactions:
    J7: -> LymphCyto ; kCyto*Cyto // Transport of Cyto to lymph Node // Max about 0.2E7
    J8: LymphCyto -> ; deltaLC*LymphCyto // Clearance of lymph cytokine
    // Amplification of T cells in Lymph Node (may have to add 0 order source)
    J9: -> LymphE ; kLE*KLLE*LymphE/(KLLE+LymphE)*LymphCyto/(KLC+LymphCyto) + aEL*(bEL-LymphE)
    J10: LymphE-> ; kLEE* LymphE //Transport of T cells to tissue
    
    // Species initializations:
    E = kLE/dE*LymphE ; // Initial number of CD8E cells in Tissue
    LymphE = aEL*bEL/(kLE+aEL) ; // Initial number of CD8E cells in Lymph Node
    Cyto = 0.0 ;
    LymphCyto = 0.0 ;
    
    // Parameter initialization:
    beta = 6.2E-5 / s_l * s_t; // virus infectivity Units, TCID^-1 /day
    k = 4.0 * s_t; // eclipse phase transition, units /day
    p = 1.0 * s_t; // virus production, units TCID/cell/day
    cV = 9.4 * s_t; // virus clearance, units /day
    deltaI2 = 1.4E-1 * s_t; // infected cell clearance, units / day
    deltaE = 5 * s_t; // infected cell clearance by CD8Es, units / day
    KdeltaE = 8.1E5 * s_v; // half saturation constant, units number of cells CD8E
    dE = 1.0 * s_t; // CD8E clearance, units /day
    cC = 1.0 * s_t; // Cytokine Production 
    deltaC = 2.0 * s_t; // Clearance of Cytokines, units of /day
    kCyto = 0.5 * s_t; // Transport of Cyto to lymph Node, units of /day
    deltaLC = 0.5 * s_t; // Clearance of lymph cytokine , units of /day
    kLE = 2.0 * s_t; // Amplification of T cells in Lymph Node (may have to add 0 order source), units of /day
    kLEE = 0.5 * s_t; // Transport of T cells to tissue, units /day
    KLC = 2E6 * s_v; // Saturation for effect of cytokine on T cell production
    KLLE = 1E7 * s_v; // Saturation of T cells in lymph node
    aEL = 1.0 * s_t;
    bEL = 4.2E5 * s_v;
    end'''
    return model_string


def immune_recruitment_model_string_ode(num_ec=1E7, num_infect=75, time_conv=1):
    model_string = f'''model {ir_model_name_ode}()
    // Scaling coefficients
    s_v = {num_ec} / 1E7;
    s_t = {time_conv};
    
    // Reactions:
    J1: TE -> I1 ; beta*TE*V //Infection rate
    J2: I1 -> I2 ; k*I1 //Eclipse Phase
    J3: I2 -> D ;  deltaI2*I2 +deltaE*E*I2/(KdeltaE+I2) // General and Immune mediated cell death
    J4: -> V; p*I2 // Virus production by virus releasing infected cells
    J4A: V -> ; cV*V // Viral clearance
    J5: -> Cyto ; cC*(I1+I2) // Cytokine Production
    J6: Cyto -> ; deltaC*Cyto // Clearance of Cytokines // Max about 0.4E7
    J7: Cyto -> LymphCyto ; kCyto*Cyto // Transport of Cyto to lymph Node // Max about 0.2E7
    J8: LymphCyto -> ; deltaLC*LymphCyto // Clearance of lymph cytokine
    // Amplification of T cells in Lymph Node
    J9: -> LymphE ; kLE*KLLE*LymphE/(KLLE+LymphE)*LymphCyto/(KLC+LymphCyto) + aEL*(bEL-LymphE)
    J10: LymphE->E ; kLEE* LymphE //Transport of T cells to tissue
    J11: E-> ; dE*E // Clearance of E cells
    
    // Species initializations:
    TE = {num_ec} ; //Intitial number of epithelal cells
    I1 = {num_infect} ; // Initial number of infected cells
    I2 = 0.0 ;
    D = 0.0 ;
    V = 0.0 ;
    E = kLE/dE*LymphE ; // Initial number of CD8E cells in Tissue
    LymphE = aEL*bEL/(kLE+aEL) ; // Initial number of CD8E cells in Lymph Node
    Cyto = 0.0 ;
    LymphCyto = 0.0 ;
    
    // Parameter initialization:
    beta = 6.2E-5 / s_v * s_t; // virus infectivity Units, TCID^-1 /day
    k = 4.0 * s_t; // eclipse phase transition, units /day
    p = 1.0 * s_t; // virus production, units TCID/cell/day
    cV = 9.4 * s_t; // virus clearance, units /day
    deltaI2 = 1.4E-1 * s_t; // infected cell clearance, units / day
    deltaE = 5 * s_t; // infected cell clearance by CD8Es, units / day
    KdeltaE = 8.1E5 * s_v; // half saturation constant, units number of cells CD8E
    dE = 1.0 * s_t; // CD8E clearance, units /day
    cC = 1.0 * s_t; // Cytokine Production 
    deltaC = 2.0 * s_t; // Clearance of Cytokines, units of /day
    kCyto = 0.5 * s_t; // Transport of Cyto to lymph Node, units of /day
    deltaLC = 0.5 * s_t; // Clearance of lymph cytokine , units of /day
    kLE = 2.0 * s_t; // Amplification of T cells in Lymph Node (may have to add 0 order source), units of /day
    kLEE = 0.5 * s_t; // Transport of T cells to tissue, units /day
    KLC = 2E6 * s_v; // Saturation for effect of cytokine on T cell production
    KLLE = 1E7 * s_v; // Saturation of T cells in lymph node
    aEL = 1.0 * s_t;
    bEL = 4.2E5 * s_v;


    // Death mechanism tracking
    -> viralDeath ; deltaI2*I2 ;
    -> cd8Death ; deltaE*E*I2/(KdeltaE+I2) ;
    viralDeath = 0;
    cd8Death = 0;

    end'''
    return model_string


def immune_recruitment_model_string_original(time_conv=1):
    model_string = f'''model {ir_model_name_orig}()
    s_t = {time_conv}

    // Reactions:
    J1: TE -> I1 ; beta*TE*V // Infection rate
    J2: I1 -> I2 ; k*I1 // Eclipse Phase
    J3: I2 -> D ;  deltaI2*I2 +deltaE*E*I2/(KdeltaE+I2) // General and Immune mediated cell death
    J4: -> V; p*I2 // Virus production by virus releasing infected cells
    J4A: V -> ; cV*V // Viral clearance
    J5: -> Cyto ; cC*(I1+I2) // Cytokine Production
    J6: Cyto -> ; deltaC*Cyto // Clearance of Cytokines // Max about 0.4E7
    J7: Cyto -> LymphCyto ; kCyto*Cyto // Transport of Cyto to lymph Node // Max about 0.2E7
    J8: LymphCyto -> ; deltaLC*LymphCyto // Clearance of lymph cytokine
    // Amplification of T cells in Lymph Node (may have to add 0 order source)
    J9: -> LymphE ; kLE*KLLE*LymphE/(KLLE+LymphE)*LymphCyto/(KLC+LymphCyto) + aEL*(bEL - LymphE)
    J10: LymphE->E ; kLEE* LymphE // Transport of T cells to tissue
    J11: E-> ; dE*E // Clearance of E cells
    
    // Species initializations:
    TE = 1E7 ; // Intitial number of epithelal cells
    I1 = 75.0 ; // Initial number of infected cells
    I2 = 0.0 ;
    D = 0.0 ;
    V = 0.0 ;
    E = kLE/dE*LymphE ; // Initial number of CD8E cells in Tissue
    LymphE = aEL*bEL/(kLE+aEL) ; // Initial number of CD8E cells in Lymph Node
    Cyto = 0.0 ;
    LymphCyto = 0.0 ;
    
    // Parameter initialization:
    beta = 6.2E-5 * s_t; // virus infectivity Units, TCID^-1 /day
    k = 4.0 * s_t; // eclipse phase transition, units /day
    p = 1.0 * s_t; // virus production, units TCID/cell/day
    cV = 9.4 * s_t; // virus clearance, units /day
    deltaI2 = 1.4E-1 * s_t; // infected cell clearance, units / day
    deltaE = 5 * s_t; // infected cell clearance by CD8Es, units / day
    KdeltaE = 8.1E5; // half saturation constant, units number of cells CD8E
    dE = 1.0 * s_t; // CD8E clearance, units /day
    cC = 1.0 * s_t; // Cytokine Production 
    deltaC = 2.0 * s_t; // Clearance of Cytokines, units of /day
    kCyto = 0.5 * s_t; // Transport of Cyto to lymph Node, units of /day
    deltaLC = 0.5 * s_t; // Clearance of lymph cytokine , units of /day
    kLE = 2.0 * s_t; // Amplification of T cells in Lymph Node (may have to add 0 order source), units of /day
    kLEE = 0.5; // Transport of T cells to tissue, units /day
    KLC = 2E6; // Saturation for effect of cytokine on T cell production
    KLLE = 1E7; // Saturation of T cells in lymph node
    aEL = 1.0 * s_t;
    bEL = 4.2E5;
    end'''
    return model_string


def ul_rate_to_prob(_rate):
    return _rate * math.exp(-_rate)
