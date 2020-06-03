import tellurium as te
import matplotlib.pyplot as plt

model_string = '''
    //Equations
    E2: -> IFN   ; P*(k11*RIGI*V+k12*V^n/(k13+V^n)+k14*IRF7P)-k21*IFN ;
    E3: -> IFNe  ; k21*IFN-t2*IFNe                                    ;
    E4: -> STATP ; P*k31*IFNe/(k32+k33*IFNe)-t3*STATP                 ;
    E5: -> IRF7  ; P*(k41*STATP+k42*IRF7P)-t4*IRF7                    ; //Why?
    E6: -> IRF7P ; P*k51*IRF7-t5*IRF7P                                ;
    E7: -> P     ; - k61*P*V                                          ; 
    E8: -> V     ; k71*V/(1+k72*IFNe)-k73*V                           ;

    //Parameters
    k11 = 0.0    ; //Validation data dependent ; zero for most cases
    k12 = 9.746  ; 
    k13 = 12.511 ; 
    k14 = 13.562 ;
    k21 = 10.385 ;
    t2  = 3.481  ;
    k31 = 45.922 ;
    k32 = 5.464  ;
    k33 = 0.068  ;
    t3  = 0.3    ;
    k41 = 0.115  ;
    k42 = 1.053  ;
    t4  = 0.3    ;
    k51 = 0.202  ;
    t5  = 0.3    ;
    k61 = 0.635  ;
    k71 = 1.537  ;
    k72 = 47.883 ;
    k73 = 0.197  ;
    n   = 3.0    ;

    //Initial Conditions
    P    = 1.0     ;   
    RIGI = 1.0     ;
    IRF7 = 0.72205 ;
    V    = 6.9e-8  ;
'''

m = te.loada(model_string)
s = m.simulate(0,36,3600)

plt.figure(figsize=(8,7))
plt.subplot(3,3,1)
plt.plot(s['time'],s['[IFN]'])
plt.title('IFN')

plt.subplot(3,3,2)
plt.plot(s['time'],s['[IFNe]'])
plt.title('IFNe')

plt.subplot(3,3,3)
plt.plot(s['time'],s['[STATP]'])
plt.title('STATP')

plt.subplot(3,3,4)
plt.plot(s['time'],s['[IRF7]'])
plt.title('IRF7')

plt.subplot(3,3,5)
plt.plot(s['time'],s['[IRF7P]'])
plt.title('IRF7P')

plt.subplot(3,3,6)
plt.plot(s['time'],s['[V]'])
plt.title('V')

plt.subplot(3,3,7)
plt.plot(s['time'],s['[P]'])
plt.ylim([0,1])
plt.title('P')

plt.tight_layout()
plt.show()