from cc3d.core.PySteppables import *
import numpy as np

plot_StandAlone = False
plot_CellModel = True
plot_SaenzModel = True

feedback = True

min_to_mcs = 10.0  # min/mcs
days_to_mcs = min_to_mcs / 1440.0  # day/mcs

'''Dynamics of Influenza Virus Infection and Pathology. 

Roberto A. Saenz, Michelle Quinlivan, Debra Elton, Shona M. Daly, Paul Digard, Ann Cullinane, Bryan T. Grenfell, John W. McCauley, James L. N. Wood and Julia R. Gog MacRae, Anthony S. Blunden, Jennifer A. Mumford, Janet

Frontiers in microbiology. J. Virol. 2010, 84(8):3974.'''

ModelString = '''        
        model Saenz()
        
//State Variables and Transitions
        E1: -> T; - beta*V*T - phi*F*T                      // T - uninfected epithelial cells
        E2: -> EP1; beta*V*T - k1*EP1                     // EP1  - 1st eclipsed phase
        E3: -> W; phi*F*T - m*beta*V*W - a*W    // W  - prerefractory state of R
        E4: -> EP2; m*beta*V*W - k2*EP2              // EP2  - IF infected, W move to this 2nd eclipsed phase 
        E5: -> R; a*W                                                // R  - refractory cells
        E6: -> I; k1*EP1 + k2*EP2 - delta*I              // I - infected cells
        E7: -> V; p*I - c*V                                         // V  - virus load
        E8: -> F; n*q*EP2 + q*I - dd*F                    // IFN released by  I cells
        E9: -> D; delta*I                                           // Dead cells
        
// Parameters from Table 1
        T0 =3.5e11       // epithelial cell population
        m = 1               // IFN-reduced infectivity
        c = 5.2             // free-virus clearance - Rate of virus clearance day^{-1}
        n = 1               // IFN-reduced production
        k1 = 2             // 1/(EP1 phase period) - days
        k2 = 2             // 1/(EP2 phase period) - days
        delta = 2        // 1/(infectious period) - days
        a = 4              // 1/(prerefractory period) - days
        
// Initial values for the average horse from Table 2
        V0av = 3.2e-1        // initial virus load - RNA copies (ml NS)^{-1}
        betaav = 1.4e-4    // infectivity rate - (RNA copies)^{-1}ml NS day^{-1}
        pav = 1.4e-5         // virus production rate - RNA copies (ml NS)^{-1} day^{-1} cell^{-}
        qav = 5e-10         // IFN production - IFN fold change day^{-1} cell^{-1}   
        phiav = 5.6e1      // IFN efficiecy - (IFN fold change)^{-1} day^{-1} 
        dav = 6.8e0         // IFN clearance - Rate of IFN clearance day^{-1}
        
// Initial Conditions
        T=T0
        beta=betaav
        phi=phiav
        EP1=0.
        W=0.
        EP2=0.
        R=0.
        I=0.
        V=V0av
        p=pav
        F=0.
        q=qav
        dd=dav
end'''

class SaenzModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Uptading max simulation steps using scaling factor to simulate 10 days
        self.get_xml_element('simulation_steps').cdata = 10.0 / days_to_mcs

        # Adding free floating antimony model
        self.add_free_floating_antimony(model_string=ModelString, model_name='Saenz',
                                        step_size=days_to_mcs)
        # Changing initial values according to discussions with Amber Smith
        # # # # state = {}
        # # # # state['I1'] = 0.0
        # # # # state['V'] = 75.0
        # # # # self.set_sbml_state(model_name='Saenz', state=state)

        # Initialize Graphic Window for Saenz ODE model
        if plot_StandAlone:
            self.plot_win = self.add_new_plot_window(title='Saenz Model Cells',
                                                 x_axis_title='Days',
                                                 y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
            self.plot_win.add_plot("T", style='Dots', color='red', size=5)
            self.plot_win.add_plot("EP1", style='Dots', color='orange', size=5)
            self.plot_win.add_plot("EP2", style='Dots', color='yellow', size=5)
            self.plot_win.add_plot("W", style='Dots', color='blue', size=5)
            self.plot_win.add_plot("R", style='Dots', color='grey', size=5)
            self.plot_win.add_plot("I", style='Dots', color='green', size=5)
            self.plot_win.add_plot("D", style='Dots', color='white', size=5)

            self.plot_win2 = self.add_new_plot_window(title='Saenz Model Virus',
                                                  x_axis_title='Days',
                                                  y_axis_title='Virus', x_scale_type='linear', y_scale_type='linear',
                                                  grid=False)
            self.plot_win2.add_plot("V", style='Dots', color='blue', size=5)
            self.plot_win2.add_plot("IFN", style='Dots', color='red', size=5)

    def step(self, mcs):
        self.timestep_sbml()
        if plot_StandAlone:
            self.plot_win.add_data_point("T", mcs * days_to_mcs,self.sbml.Saenz['T'] / self.sbml.Saenz['T0'])
            self.plot_win.add_data_point("EP1", mcs * days_to_mcs,self.sbml.Saenz['EP1'] / self.sbml.Saenz['T0'])
            self.plot_win.add_data_point("EP2", mcs * days_to_mcs,self.sbml.Saenz['EP2'] / self.sbml.Saenz['T0'])
            self.plot_win.add_data_point("W", mcs * days_to_mcs,self.sbml.Saenz['W'] / self.sbml.Saenz['T0'])
            self.plot_win.add_data_point("R", mcs * days_to_mcs,self.sbml.Saenz['R'] / self.sbml.Saenz['T0'])
            self.plot_win.add_data_point("I", mcs * days_to_mcs,self.sbml.Saenz['I'] / self.sbml.Saenz['T0'])
            self.plot_win2.add_data_point("V", mcs * days_to_mcs, np.log10(self.sbml.Saenz['V']))
            if  self.sbml.Saenz['F'] > .1:
                self.plot_win2.add_data_point("IFN", mcs * days_to_mcs, np.log10(self.sbml.Saenz['F']))
            self.plot_win2.add_data_point("D", mcs * days_to_mcs, np.log10(self.sbml.Saenz['D']))


class CellularModelSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # set initial model parameters
        self.initial_uninfected = len(self.cell_list)  # Scale factor for fraction of cells infected
        self.ExtracellularVirus = self.sbml.Saenz['V']
        self.IFNReleased = self.sbml.Saenz['F']
        # Initial number of I cells
        Iinit=50
        count=Iinit
        while count > 0:
                xrand=np.random.random()*self.dim.x
                yrand=np.random.random()*self.dim.y
                cell = self.cell_field[xrand,yrand, 0]
                if cell:
                    if cell.type==self.U:
                        cell.type = self.EP1
                        count-=1
        self.sbml.Saenz['EP1']=Iinit*self.sbml.Saenz['T0']/len(self.cell_list)
        # Initial number of W cells
        Winit=50
        count=Winit
        while count > 0:
                xrand=np.random.random()*self.dim.x
                yrand=np.random.random()*self.dim.y
                cell = self.cell_field[xrand,yrand, 0]
                if cell:
                       if cell.type==self.U:
                                cell.type = self.W
                                count-=1
        self.sbml.Saenz['W']=Winit*self.sbml.Saenz['T0']/len(self.cell_list)
        #print('----------------------------------------------------------------------IniVal = ',self.ExtracellularVirus,self.IFNReleased)

        self.plot_win3 = self.add_new_plot_window(title='CPM Cells',
                                                  x_axis_title='days',
                                                  y_axis_title='Variables', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)
        self.plot_win3.add_plot("U", style='Lines', color='red', size=5)
        self.plot_win3.add_plot("EP1", style='Lines', color='orange', size=5)
        self.plot_win3.add_plot("EP2", style='Lines', color='yellow', size=5)
        self.plot_win3.add_plot("W", style='Lines', color='blue', size=5)
        self.plot_win3.add_plot("R", style='Lines', color='grey', size=5)
        self.plot_win3.add_plot("I", style='Lines', color='green', size=5)
        self.plot_win3.add_plot("D", style='Lines', color='white', size=5)

        self.plot_win3.add_plot("AU", style='Dots', color='red', size=5)
        self.plot_win3.add_plot("AEP1", style='Dots', color='orange', size=5)
        self.plot_win3.add_plot("AEP2", style='Dots', color='yellow', size=5)
        self.plot_win3.add_plot("AW", style='Dots', color='blue', size=5)
        self.plot_win3.add_plot("AR", style='Dots', color='grey', size=5)
        self.plot_win3.add_plot("AI", style='Dots', color='green', size=5)
        self.plot_win3.add_plot("AD", style='Dots', color='white', size=5)

        self.plot_win4 = self.add_new_plot_window(title='CPM Virus',
                                                  x_axis_title='days',
                                                  y_axis_title='Variables', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)
        self.plot_win4.add_plot("V", style='Lines', color='blue', size=5)
        self.plot_win4.add_plot("F", style='Lines', color='red', size=5)
        self.plot_win4.add_plot("AV", style='Dots', color='blue', size=5)
        self.plot_win4.add_plot("AF", style='Dots', color='red', size=5)

    def step(self, mcs):
        # Transition rule from U to EP1
        beta = self.sbml.Saenz['beta'] * days_to_mcs
        fi = self.sbml.Saenz['phi'] * days_to_mcs
        if feedback ==  True:
            V = self.ExtracellularVirus
            F = self.IFNReleased
        else:
            V = self.sbml.Saenz['V']
            F = self.sbml.Saenz['F']
        
        # 1st Transition rule for U: -> EP1 cells
        p_UtoEP1 = beta * V
        print('---------',beta,V,p_UtoEP1)
        for cell in self.cell_list_by_type(self.U):
            if np.random.random() < (p_UtoEP1 ): 
                cell.type = self.EP1

        # 2nd Transition rule for U: -> W cells
        p_UtoW = fi*F
        print(fi,F,p_UtoW)
        for cell in self.cell_list_by_type(self.U):
            if np.random.random() < (p_UtoW ): 
                cell.type = self.W

        # Transition rule for EP1: -> I cells
        k1 = self.sbml.Saenz['k1'] * days_to_mcs
        p_EP1toI = k1
        for cell in self.cell_list_by_type(self.EP1):
            if np.random.random() < p_EP1toI:
                cell.type = self.I

        # 1st Transition rule for  W : -> EP2
        m = self.sbml.Saenz['m']
        p_WtoEP2 = m*beta*V
        for cell in self.cell_list_by_type(self.W):
            if np.random.random() < ( p_WtoEP2 ): 
                cell.type = self.EP2

        # 2nd Transition rule for  W to R
        a = self.sbml.Saenz['a'] * days_to_mcs
        p_WtoR = a
        for cell in self.cell_list_by_type(self.W):
            if np.random.random() < ( p_WtoR ): 
                cell.type = self.R

        # 2nd Transition rule for I: EP2 cells
        k2 = self.sbml.Saenz['k2'] * days_to_mcs
        p_EP2toI = k2
        for cell in self.cell_list_by_type(self.EP2):
            if np.random.random() < p_EP2toI:
                cell.type = self.I
 
        # Transition rule for I                         
        delta= self.sbml.Saenz['delta'] * days_to_mcs
        p_ItodeltaI = delta
        for cell in self.cell_list_by_type(self.I):
            if np.random.random() < p_ItodeltaI :
                cell.type = self.D

        # Extracellular Virus
        V = self.ExtracellularVirus
        p = self.sbml.Saenz['p'] / self.initial_uninfected * self.sbml.Saenz['T0'] * days_to_mcs
        c = self.sbml.Saenz['c'] * days_to_mcs
        for cell in self.cell_list_by_type(self.I):
            self.ExtracellularVirus += p
        self.ExtracellularVirus -= c * V

        # Extracellular IFN
        # F = self.IFNReleased
        q = self.sbml.Saenz['q'] / self.initial_uninfected * self.sbml.Saenz['T0'] * days_to_mcs
        d = self.sbml.Saenz['dd'] * days_to_mcs
        n = self.sbml.Saenz['n']
        for cell in self.cell_list_by_type(self.EP2):
            self.IFNReleased += n*q  # ?????
        for cell in self.cell_list_by_type(self.I):
            self.IFNReleased += q       # ?????
        self.IFNReleased -=d*self.IFNReleased

        if plot_CellModel:
            #print('*******',self.ExtracellularVirus,self.IFNReleased)
            self.plot_win3.add_data_point("U", mcs * days_to_mcs, len(self.cell_list_by_type(self.U)) / self.initial_uninfected)
            self.plot_win3.add_data_point("EP1", mcs * days_to_mcs,len(self.cell_list_by_type(self.EP1)) / self.initial_uninfected)
            self.plot_win3.add_data_point("EP2", mcs * days_to_mcs,len(self.cell_list_by_type(self.EP2)) / self.initial_uninfected)
            self.plot_win3.add_data_point("W", mcs * days_to_mcs, len(self.cell_list_by_type(self.W)) / self.initial_uninfected)
            self.plot_win3.add_data_point("R", mcs * days_to_mcs,len(self.cell_list_by_type(self.R)) / self.initial_uninfected)
            self.plot_win3.add_data_point("I", mcs * days_to_mcs,len(self.cell_list_by_type(self.I)) / self.initial_uninfected)
            self.plot_win3.add_data_point("D", mcs * days_to_mcs,len(self.cell_list_by_type(self.D)) / self.initial_uninfected)
            self.plot_win4.add_data_point("V", mcs * days_to_mcs, np.log10(self.ExtracellularVirus))
            if self.IFNReleased > 0.1:
                print('VtoF_CM =',self.ExtracellularVirus/self.IFNReleased)
                self.plot_win4.add_data_point("F", mcs * days_to_mcs, np.log10(F))                  #self.IFNReleased))

        if plot_SaenzModel:
            print('VtoF =',self.sbml.Saenz['V']/self.sbml.Saenz['F'])
            self.plot_win3.add_data_point("AU", mcs * days_to_mcs, self.sbml.Saenz['T'] / self.sbml.Saenz['T0'])
            self.plot_win3.add_data_point("AEP1", mcs * days_to_mcs,self.sbml.Saenz['EP1'] / self.sbml.Saenz['T0'])
            self.plot_win3.add_data_point("AEP2", mcs * days_to_mcs,self.sbml.Saenz['EP2'] / self.sbml.Saenz['T0'])
            self.plot_win3.add_data_point("AW", mcs * days_to_mcs, self.sbml.Saenz['W'] / self.sbml.Saenz['T0'])
            self.plot_win3.add_data_point("AR", mcs * days_to_mcs,self.sbml.Saenz['R'] / self.sbml.Saenz['T0'])
            self.plot_win3.add_data_point("AI", mcs * days_to_mcs,self.sbml.Saenz['I'] / self.sbml.Saenz['T0'])
            self.plot_win3.add_data_point("AD", mcs * days_to_mcs,self.sbml.Saenz['D'] / self.sbml.Saenz['T0'])
            self.plot_win4.add_data_point("AV", mcs * days_to_mcs, np.log10(self.sbml.Saenz['V']))
            if  self.sbml.Saenz['F'] > .1:
                self.plot_win4.add_data_point("AF", mcs * days_to_mcs, np.log10(self.sbml.Saenz['F']))


class StatisticsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.initial_uninfected = len(self.cell_list)

        self.cellular_infection =  False
        self.cellular_infection_time = 0.0
        self.Saenzmodel_infection = False
        self.Saenzmodel_infection_time = 0.0
        self.infection_threshold = 0.1

        self.plot_win5 = self.add_new_plot_window(title='Residuals',
                                                  x_axis_title='days',
                                                  y_axis_title='Variables', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)
        self.plot_win5.add_plot("dU", style='Lines', color='red', size=5)
        self.plot_win5.add_plot("dEP1", style='Lines', color='orange', size=5)
        self.plot_win5.add_plot("dEP2", style='Lines', color='yellow', size=5)
        self.plot_win5.add_plot("dW", style='Lines', color='blue', size=5)
        self.plot_win5.add_plot("dR", style='Lines', color='grey', size=5)
        self.plot_win5.add_plot("dI", style='Lines', color='green', size=5)
        self.plot_win5.add_plot("dD", style='Lines', color='white', size=5)

    def step(self, mcs):
        if self.cellular_infection == False:
            if len(self.cell_list_by_type(self.EP1))/self.initial_uninfected >= self.infection_threshold:
                self.cellular_infection_time = mcs
                self.cellular_infection = True

        if self.Saenzmodel_infection == False:
            if self.sbml.Saenz['EP1']/self.sbml.Saenz['T0'] >= self.infection_threshold:
                self.Saenzmodel_infection_time = mcs
                self.Saenzmodel_infection =  True

        print("Cellular Infection = ", self.cellular_infection_time * days_to_mcs)
        print("ODE Infection = ", self.Saenzmodel_infection_time * days_to_mcs)

        dU = (len(self.cell_list_by_type(self.U)) / self.initial_uninfected) - (self.sbml.Saenz['T'] / self.sbml.Saenz['T0'])
        dEP1 = (len(self.cell_list_by_type(self.EP1)) / self.initial_uninfected) - (self.sbml.Saenz['EP1'] / self.sbml.Saenz['T0'])
        dEP2 = (len(self.cell_list_by_type(self.EP2)) / self.initial_uninfected) - (self.sbml.Saenz['EP2'] / self.sbml.Saenz['T0'])
        dW = (len(self.cell_list_by_type(self.W)) / self.initial_uninfected) - (self.sbml.Saenz['W'] / self.sbml.Saenz['T0'])
        dR = (len(self.cell_list_by_type(self.R)) / self.initial_uninfected) - (self.sbml.Saenz['R'] / self.sbml.Saenz['T0'])
        dI = (len(self.cell_list_by_type(self.I)) / self.initial_uninfected) - (self.sbml.Saenz['I'] / self.sbml.Saenz['T0'])
        dD = (len(self.cell_list_by_type(self.D)) / self.initial_uninfected) - (self.sbml.Saenz['D'] / self.sbml.Saenz['T0'])

        self.plot_win5.add_data_point("dU", mcs * days_to_mcs, dU)
        self.plot_win5.add_data_point("dEP1", mcs * days_to_mcs, dEP1)
        self.plot_win5.add_data_point("dEP2", mcs * days_to_mcs, dEP2)
        self.plot_win5.add_data_point("dW", mcs * days_to_mcs, dW)
        self.plot_win5.add_data_point("dR", mcs * days_to_mcs, dR)
        self.plot_win5.add_data_point("dI", mcs * days_to_mcs, dI)
        self.plot_win5.add_data_point("dD", mcs * days_to_mcs, dD)

#         # Plot lagged differences between cell populations
#         # Start when both populations are infected
#         # Josh--you will need to save the time series in a set of lists so you can do this
#         # Once both times series have had infection begin
#         # Plot (x(t-starttime)-X(t-cellstarttime)) for each series
#         # Could do same thing to show virus with lags and also RMS deviation with lags
#         #Next step, have virus diffuse and cells infected by viral field rather than the external variable
