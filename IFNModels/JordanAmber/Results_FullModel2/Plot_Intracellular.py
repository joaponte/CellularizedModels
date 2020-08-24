import numpy as np
import matplotlib.pyplot as plt

CC3Dnames = ['Time','CC3DV','CC3DH','CC3DP','CC3DIFNe_Scalar','CC3DIFNe_Field','CC3DSTATP','CC3DIRF7','CC3DIRF7P','CC3DIFN']
CC3Dts = 262

#Plot Virus
V = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    V[:,r-1] = f['CC3DV']
    plt.plot(T, f['CC3DV'],color='#6600FF',linewidth = 2.0, alpha=0.1)
plt.plot(T,np.mean(V,1),color='#6600FF',linewidth = 3.0,label='Virus')
plt.title('Virus Level')
plt.xlabel('Time (hrs)')
plt.ylabel('Virus Level')
plt.legend(loc=1)
plt.savefig('Fig.V.pdf')
plt.clf()

# Plot Cell Viability
H = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    H[:,r-1] = f['CC3DH']
    plt.plot(T, f['CC3DH'],color='#339900',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(H,1),color='#339900',linewidth = 3.0,label='Cell Viability')
plt.title('Cell Viability')
plt.xlabel('Time (hrs)')
plt.ylabel('Viability/cell')
plt.legend(loc=1)
plt.savefig('Fig.H.pdf')
plt.clf()

# Plot Live Cells
P = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    P[:,r-1] = f['CC3DP']
    plt.plot(T, f['CC3DP'],color='red',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(P,1),color='red',linewidth = 3.0,label='Live Cells')
plt.title('Fraction of Live Cells')
plt.xlabel('Time (hrs)')
plt.ylabel('Fraction of Live Cells')
plt.legend(loc=1)
plt.savefig('Fig.P.pdf')
plt.clf()

# Plot IFNe
IFNe = np.zeros((CC3Dts, 20))
IFNe_Scalar = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IFNe[:,r-1] = f['CC3DIFNe_Field']
    IFNe_Scalar[:,r-1] = f['CC3DIFNe_Scalar']
    plt.plot(T, f['CC3DIFNe_Field'],color='#993300',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IFNe,1),color='#993300',linewidth = 3.0,label='IFNe')
plt.title('Extracellular IFN')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IFNe] $\mu$M/cell')
plt.ylim([0,1.15*max(np.mean(IFNe,1))])
plt.legend(loc=1)
plt.savefig('Fig.IFNe.pdf')
plt.clf()

# Plot STATP
STATP = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    STATP[:,r-1] = f['CC3DSTATP']
    plt.plot(T, f['CC3DSTATP'],color='#0033FF',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(STATP,1),color='#0033FF',linewidth = 3.0,label='STATP')
plt.title('STATP')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[STATP] $\mu$M/cell')
plt.ylim([0,1.15*max(np.mean(STATP,1))])
plt.legend(loc=1)
plt.savefig('Fig.STATP.pdf')
plt.clf()

# Plot IRF7
IRF7 = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IRF7[:,r-1] = f['CC3DIRF7']
    plt.plot(T, f['CC3DIRF7'],color='#FF6600',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IRF7,1),color='#FF6600',linewidth = 3.0,label='IRF7')
plt.title('IRF7')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IRF7] $\mu$M/cell')
plt.ylim([0,1.15*max(np.mean(IRF7,1))])
plt.legend(loc=1)
plt.savefig('Fig.IRF7.pdf')
plt.clf()

# Plot IRF7P
IRF7P = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IRF7P[:,r-1] = f['CC3DIRF7P']
    plt.plot(T, f['CC3DIRF7P'],color='#666666',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IRF7P,1),color='#666666',linewidth = 3.0,label='IRF7P')
plt.title('IRF7P')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IRF7P] $\mu$M/cell')
plt.ylim([0,1.15*max(np.mean(IRF7P,1))])
plt.legend(loc=1)
plt.savefig('Fig.IRF7P.pdf')
plt.clf()

# Plot INF
IFN = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IFN[:,r-1] = f['CC3DIFN']
    plt.plot(T, f['CC3DIFN'],color='#CC0033',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IFN,1),color='#CC0033',linewidth = 3.0,label='IFN')
plt.title('IFN')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IFN] $\mu$M/cell')
plt.ylim([0,1.15*max(np.mean(IFN,1))])
plt.legend(loc=1)
plt.savefig('Fig.IFN.pdf')
plt.clf()