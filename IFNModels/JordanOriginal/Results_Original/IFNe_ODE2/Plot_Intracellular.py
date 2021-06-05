import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

CC3Dnames = ['Time','CC3DV','CC3DH','CC3DP','CC3DIFNe_Scalar','CC3DIFNe_Field','CC3DSTATP','CC3DIRF7','CC3DIRF7P','CC3DIFN']
CC3Dts = 180

ODEnames = ['Time','ODEV','ODEH','ODEP','ODEIFNe','ODESTATP','ODEIRF7','ODEIRF7P','ODEIFN']
ODEf = np.genfromtxt('JordanOriginalODE_1.txt', skip_header=1, delimiter=',', names=ODEnames)

#Plot Virus
V = np.zeros((CC3Dts, 20))
plt.figure(figsize=([6.2,4.8]))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    V[:,r-1] = f['CC3DV']
    plt.plot(T, f['CC3DV'],color='#6600FF',linewidth = 4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODEV'][::2],'.',color='#6600FF',markersize=10,label='ODE')
plt.plot(T,np.mean(V,1),color='#6600FF',linewidth = 4.0,label='ABM')
plt.title('Virus Level')
plt.xlabel('Time (hrs)')
plt.ylabel('Virus Level')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.V.pdf',transparent=True)
plt.clf()

# Plot Cell Viability
H = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    H[:,r-1] = f['CC3DH']
    plt.plot(T, f['CC3DH'],color='#339900',linewidth = 4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODEH'][::2],'.',color='#339900',markersize=10,label='ODE')
plt.plot(T,np.mean(H,1),color='#339900',linewidth = 4.0,label='ABM')
plt.title('Cell Viability')
plt.xlabel('Time (hrs)')
plt.ylabel('Viability/cell')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.H.pdf',transparent=True)
plt.clf()

# Plot Live Cells
P = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    P[:,r-1] = f['CC3DP']
    plt.plot(T, f['CC3DP'],color='red',linewidth = 4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODEP'][::2],'.',color='red',markersize=10,label='ODE')
plt.plot(T,np.mean(P,1),color='red',linewidth = 4.0,label='ABM')
plt.title('Fraction of Live Cells')
plt.xlabel('Time (hrs)')
plt.ylabel('Fraction of Live Cells')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.P.pdf',transparent=True)
plt.clf()

# Plot IFNe
IFNe = np.zeros((CC3Dts, 20))
IFNe_Scalar = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IFNe[:,r-1] = f['CC3DIFNe_Field']
    IFNe_Scalar[:,r-1] = f['CC3DIFNe_Scalar']
    plt.plot(T, f['CC3DIFNe_Field'],color='#993300',linewidth = 4.0, alpha=0.15)
    plt.plot(T, f['CC3DIFNe_Scalar'],'--', color='#993300', linewidth=4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODEIFNe'][::2],'.',color='#993300',markersize=10,label='ODE')
plt.plot(T,np.mean(IFNe,1),color='#993300',linewidth = 4.0,label='ABM_Field')
plt.plot(T,np.mean(IFNe_Scalar,1),'--', color='#993300',linewidth = 4.0,label='ABM_Scalar')
plt.title('Extracellular IFN')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IFNe] $\mu$M/cell')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.IFNe.pdf',transparent=True)
plt.clf()

# Plot STATP
STATP = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    STATP[:,r-1] = f['CC3DSTATP']
    plt.plot(T, f['CC3DSTATP'],color='#0033FF',linewidth = 4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODESTATP'][::2],'.',color='#0033FF',markersize=10,label='ODE')
plt.plot(T,np.mean(STATP,1),color='#0033FF',linewidth = 4.0,label='ABM')
plt.title('STATP')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[STATP] $\mu$M/cell')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.STATP.pdf',transparent=True)
plt.clf()

# Plot IRF7
IRF7 = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IRF7[:,r-1] = f['CC3DIRF7']
    plt.plot(T, f['CC3DIRF7'],color='#FF6600',linewidth = 4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODEIRF7'][::2],'.',color='#FF6600',markersize=10,label='ODE')
plt.plot(T,np.mean(IRF7,1),color='#FF6600',linewidth = 4.0,label='ABM')
plt.title('IRF7')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IRF7] $\mu$M/cell')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.IRF7.pdf',transparent=True)
plt.clf()

# Plot IRF7P
IRF7P = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IRF7P[:,r-1] = f['CC3DIRF7P']
    plt.plot(T, f['CC3DIRF7P'],color='#666666',linewidth = 4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODEIRF7P'][::2],'.',color='#666666',markersize=10,label='ODE')
plt.plot(T,np.mean(IRF7P,1),color='#666666',linewidth = 4.0,label='ABM')
plt.title('IRF7P')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IRF7P] $\mu$M/cell')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.IRF7P.pdf',transparent=True)
plt.clf()

# Plot INF
IFN = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('JordanOriginalCC3D_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IFN[:,r-1] = f['CC3DIFN']
    plt.plot(T, f['CC3DIFN'],color='#CC0033',linewidth = 4.0, alpha=0.15)
plt.plot(ODEf['Time'][::2], ODEf['ODEIFN'][::2],'.',color='#CC0033',markersize=10,label='ODE')
plt.plot(T,np.mean(IFN,1),color='#CC0033',linewidth = 4.0,label='ABM')
plt.title('IFN')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IFN] $\mu$M/cell')
plt.legend(loc=1)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.IFN.pdf',transparent=True)
plt.clf()