import numpy as np
import matplotlib.pyplot as plt

CellNames = ['Time','U','I1','I2','D','Ve']
ts = 300
replicates = 20
U = np.zeros((ts, replicates))
I1 = np.zeros((ts, replicates))
I2 = np.zeros((ts, replicates))
D = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    plt.plot(f['Time'], f['U'], color='blue',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['I1'], color='orange',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['I2'], color='red',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['D'], color='purple',linewidth = 3.0, alpha=0.075)
    U[:, r - 1] = f['U']
    I1[:, r - 1] = f['I1']
    I2[:, r - 1] = f['I2']
    D[:, r - 1] = f['D']
plt.plot(f['Time'],np.mean(U,1),color='blue',linewidth = 3.0,label='U')
plt.plot(f['Time'],np.mean(I1,1),color='orange',linewidth = 3.0,label='I1')
plt.plot(f['Time'],np.mean(I2,1),color='red',linewidth = 3.0,label='I2')
plt.plot(f['Time'],np.mean(D,1),color='purple',linewidth = 3.0,label='D')
plt.legend(loc=1)
plt.savefig('Fig.Cells.pdf')
plt.clf()

Ve = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    plt.plot(f['Time'], f['Ve'], color='black',linewidth = 4.0, alpha=0.1)
    Ve[:, r - 1] = f['Ve']
plt.plot(f['Time'],np.mean(Ve,1),color='black',linewidth = 3.0,label='Viral Load')
plt.ylim([1,None])
plt.legend(loc=2)
plt.yscale('log')
plt.savefig('Fig.ViralLoad.pdf')
