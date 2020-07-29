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

    //Equations
    E8: -> Cl ; kc * C // Cytokine transport to the lymph node
    E9: Cl -> ; ccl * Cl // Lymph node cytokine decay
    E10: -> El ; pel*Cl + rel*Kel*Cl*El/(Kel + El) ;
    E11: El -> ; ke*El   ;
    
    
    //Parameters
    beta = 6.2E-5 * s_t / s_l ; // 1.0/(TCID*day)
    k = 4.0 * s_t ; // 1.0/day
    p = 1.0 * s_t ; // TCID/cell/day
    c = 9.4 * s_t ; // 1.0/day
    d = 2.4E-1 * s_t; // 1.0/day
    pc = 1.0 * s_t ;
    cc = 2.0 * s_t ; // 1.0/day
    kc = 0.5 * s_t ; // 1.0/day
    ccl = 0.5 * s_t ; // 1.0/day
    pel = 1E-4 * s_t ;
    rel = 1.0E-7 * s_t / s_v ;
    Kel = 1E3 * s_v ;
    ke = 0.5 * s_t ;
    dE = 1.0 * s_t ;
    dei2 = 12E-2 * s_t / s_v ; //15E3
    kei2 = 5E4 * s_v ;
    
    //Inputs
    V = 0.0
    
    //Initial Conditions
    T0 = {num_ec}
    T = T0
    end'''
    return model_string


def immune_recruitment_model_string_ode(num_ec=1E7, num_infect=75, time_conv=1):
    model_string = f'''model {ir_model_name_ode}()
    // Scaling coefficients
    s_v = {num_ec} / 1E7;
    s_t = {time_conv};
    
    //Equations
    E1: T -> I1 ; beta*T*V ; //Infection rate
    E2: I1 -> I2 ; k*I1 ; //Infection rate
    E3: -> V ; p*I2 ; //Virus Production
    E4: V -> ; c*V ; // Virus Decay
    E5: I2 -> D ; d*I2 ; // Infected Cell Clearance (apopotosis + innate immune response)
    E6: -> C ; pc*(I1+I2) ; // Cytokine production
    E7: C -> ; cc*C ; // Cytokine decay
    E8: C -> Cl ; kc * C // Cytokine transport to the lymph node
    E9: Cl -> ; ccl * Cl // Lymph node cytokine decay
    E10: -> El ; pel*Cl + rel*Kel*Cl*El/(Kel + El) ;
    E11: El -> E; ke*El   ;
    E12: E ->  ; dE*E     ;
    E13: I2 -> D ; kei2*dei2*E/(E+kei2+I2)*I2 ;
    
    
    //Parameters
    beta = 6.2E-5 * s_t / s_v ; // 1.0/(TCID*day)
    k = 4.0 * s_t ; // 1.0/day
    p = 1.0 * s_t ; // TCID/cell/day
    c = 9.4 * s_t ; // 1.0/day
    d = 2.4E-1 * s_t; // 1.0/day
    pc = 1.0 * s_t ;
    cc = 2.0 * s_t ; // 1.0/day
    kc = 0.5 * s_t ; // 1.0/day
    ccl = 0.5 * s_t ; // 1.0/day
    pel = 1E-4 * s_t ;
    rel = 1.0E-7 * s_t / s_v ;
    Kel = 1E3 * s_v ;
    ke = 0.5 * s_t ;
    dE = 1.0 * s_t ;
    dei2 = 12E-2 * s_t / s_v ; //15E3
    kei2 = 5E4 * s_v ;
    
    //Inputs
    V = 0.0
    
    //Initial Conditions
    T0 = {num_ec}
    T = T0
    I1 = {num_infect}


    // Death mechanism tracking
    -> viralDeath ; d*I2 ;
    -> cd8Death ; kei2*dei2*E/(E+kei2+I2)*I2 ;
    viralDeath = 0;
    cd8Death = 0;

    end'''
    return model_string


def immune_recruitment_model_string_original(time_conv=1):
    model_string = f'''model {ir_model_name_orig}()
    s_t = {time_conv}

    //Equations
    E1: T -> I1 ; beta*T*V ; //Infection rate
    E2: I1 -> I2 ; k*I1 ; //Infection rate
    E3: -> V ; p*I2 ; //Virus Production
    E4: V -> ; c*V ; // Virus Decay
    E5: I2 -> D ; d*I2 ; // Infected Cell Clearance (apopotosis + innate immune response)
    E6: -> C ; pc*(I1+I2) ; // Cytokine production
    E7: C -> ; cc*C ; // Cytokine decay
    E8: C -> Cl ; kc * C // Cytokine transport to the lymph node
    E9: Cl -> ; ccl * Cl // Lymph node cytokine decay
    E10: -> El ; pel*Cl + rel*Kel*Cl*El/(Kel + El) ;
    E11: El -> E; ke*El   ;
    E12: E ->  ; dE*E     ;
    E13: I2 -> D ; kei2*dei2*E/(E+kei2+I2)*I2 ;
    
    
    //Parameters
    beta = 6.2E-5 * s_t ; // 1.0/(TCID*day)
    k = 4.0 * s_t ; // 1.0/day
    p = 1.0 * s_t ; // TCID/cell/day
    c = 9.4 * s_t ; // 1.0/day
    d = 2.4E-1 * s_t; // 1.0/day
    pc = 1.0 * s_t ;
    cc = 2.0 * s_t ; // 1.0/day
    kc = 0.5 * s_t ; // 1.0/day
    ccl = 0.5 * s_t ; // 1.0/day
    pel = 1E-4 * s_t ;
    rel = 1.0E-7 * s_t ;
    Kel = 1E3 ;
    ke = 0.5 * s_t ;
    dE = 1.0 * s_t ;
    dei2 = 12E-2 * s_t ; //15E3
    kei2 = 5E4 ;
    
    //Inputs
    V = 0.0
    
    //Initial Conditions
    T0 = 1E7
    T = T0
    I1 = 75 
    end'''
    return model_string


def ul_rate_to_prob(_rate):
    return _rate * math.exp(-_rate)
