import numpy as np
import matplotlib.pyplot as plt

CellNames = ['Time','U','I1','I2','D','Ve','IFNe']
ts = 262
replicates = 20
U = np.zeros((ts, replicates))
I1 = np.zeros((ts, replicates))
I2 = np.zeros((ts, replicates))
D = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    #plt.plot(f['Time'], f['U'], color='blue',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['I1'], color='orange',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['I2'], color='red',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['D'], color='purple',linewidth = 3.0, alpha=0.075)
    U[:, r - 1] = f['U']
    I1[:, r - 1] = f['I1']
    I2[:, r - 1] = f['I2']
    D[:, r - 1] = f['D']
#plt.plot(f['Time'],np.mean(U,1),color='blue',linewidth = 3.0,label='U')
plt.plot(f['Time'],np.mean(I1,1),color='orange',linewidth = 3.0,label='I1')
plt.plot(f['Time'],np.mean(I2,1),color='red',linewidth = 3.0,label='I2')
plt.plot(f['Time'],np.mean(D,1),color='purple',linewidth = 3.0,label='D')
plt.title('Cells Types')
plt.xlabel('Time (hrs)')
plt.ylabel('Fraction of Cells by Type')
plt.legend(loc=2)
plt.tight_layout()
plt.savefig('Fig.Cells.pdf')
plt.clf()

Ve = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    plt.plot(f['Time'], f['Ve'], color='#666699',linewidth = 4.0, alpha=0.1)
    Ve[:, r - 1] = f['Ve']
plt.plot(f['Time'],np.mean(Ve,1),color='#666699',linewidth = 3.0,label='Viral Load')
plt.title('Viral Load')
plt.xlabel('Time (hrs)')
plt.ylabel(r'log$_{10}$ PFU/ml')
plt.yscale('log')
plt.ylim([1,None])
plt.legend(loc=2)
plt.tight_layout()
plt.savefig('Fig.ViralLoad.pdf')
plt.clf()

IFNe = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    plt.plot(f['Time'], f['IFNe'], color='#993300',linewidth = 4.0, alpha=0.1)
    IFNe[:, r - 1] = f['IFNe']
plt.plot(f['Time'],np.mean(IFNe,1),color='#993300',linewidth = 3.0,label='IFNe')
plt.title('Extracellular IFN')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IFNe] $\mu$M')
plt.yscale('log')
plt.legend(loc=2)
plt.tight_layout()
plt.savefig('Fig.TotalIFNe.pdf')
