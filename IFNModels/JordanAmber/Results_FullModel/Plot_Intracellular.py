import numpy as np
import matplotlib.pyplot as plt

CC3Dnames = ['Time','CC3DV','CC3DH','CC3DP','CC3DIFNe_Scalar','CC3DIFNe_Field','CC3DSTATP','CC3DIRF7','CC3DIRF7P','CC3DIFN']
CC3Dts = 300

#Plot Virus
V = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    V[:,r-1] = f['CC3DV']
    plt.plot(T, f['CC3DV'],color='purple',linewidth = 2.0, alpha=0.1)
plt.plot(T,np.mean(V,1),color='purple',linewidth = 3.0,label='ABM')
plt.title('Virus Level')
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
    plt.plot(T, f['CC3DH'],color='green',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(H,1),color='green',linewidth = 3.0,label='ABM')
plt.title('Cell Viability')
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
plt.plot(T,np.mean(P,1),color='red',linewidth = 3.0,label='ABM')
plt.title('Fraction of Cells Alive')
plt.legend(loc=1)
plt.savefig('Fig.P.pdf')
plt.clf()

# Plot Live Cells
IFNe = np.zeros((CC3Dts, 20))
for r in range(1,21):
    f = np.genfromtxt('FullModelIntracellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CC3Dnames,
                      max_rows = CC3Dts)
    T = f['Time']
    IFNe[:,r-1] = f['CC3DIFNe_Field']
    plt.plot(T, f['CC3DIFNe_Field'],color='orange',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IFNe,1),color='orange',linewidth = 3.0,label='ABM')
plt.title('Extracellular IFN')
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
    plt.plot(T, f['CC3DSTATP'],color='blue',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(STATP,1),color='blue',linewidth = 3.0,label='ABM')
plt.title('STATP')
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
    plt.plot(T, f['CC3DIRF7'],color='magenta',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IRF7,1),color='magenta',linewidth = 3.0,label='ABM')
plt.title('IRF7')
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
    plt.plot(T, f['CC3DIRF7P'],color='brown',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IRF7P,1),color='brown',linewidth = 3.0,label='ABM')
plt.title('IRF7P')
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
    plt.plot(T, f['CC3DIFN'],color='black',linewidth = 2.0, alpha=0.15)
plt.plot(T,np.mean(IFN,1),color='black',linewidth = 3.0,label='ABM')
plt.title('IFN')
plt.legend(loc=1)
plt.savefig('Fig.IFN.pdf')
plt.clf()