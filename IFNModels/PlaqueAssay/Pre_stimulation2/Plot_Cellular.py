import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

CellNames = ['Time','U','I1','I2','D','Ve','IFNe']
ts = 480
replicates = 20
U = np.zeros((ts, replicates))
I1 = np.zeros((ts, replicates))
I2 = np.zeros((ts, replicates))
D = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    #plt.plot(f['Time'], f['U'], color='blue',linewidth = 3.0, alpha=0.075)
    # plt.plot(f['Time'], f['I1'], color='orange',linewidth = 3.0, alpha=0.075)
    # plt.plot(f['Time'], f['I2'], color='red',linewidth = 3.0, alpha=0.075)
    # plt.plot(f['Time'], f['D'], color='purple',linewidth = 3.0, alpha=0.075)
    U[:, r - 1] = f['U']
    I1[:, r - 1] = f['I1']
    I2[:, r - 1] = f['I2']
    D[:, r - 1] = f['D']
# plt.plot(f['Time'],np.mean(U,1),color='blue',linewidth = 3.0,label='U')
plt.plot(f['Time'],np.mean(I1,1),color='orange',linewidth = 3.0,label='I1')
plt.fill_between(f['Time'], np.amin(I1, 1), np.amax(I1, 1), facecolor='orange', interpolate=True, alpha=0.2)
plt.plot(f['Time'],np.mean(I2,1),color='red',linewidth = 3.0,label='I2')
plt.fill_between(f['Time'], np.amin(I2, 1), np.amax(I2, 1), facecolor='red', interpolate=True, alpha=0.2)
plt.plot(f['Time'],np.mean(D,1),color='purple',linewidth = 3.0,label='D')
plt.fill_between(f['Time'], np.amin(D, 1), np.amax(D, 1), facecolor='purple', interpolate=True, alpha=0.2)
plt.title('Cells Types')
plt.xlabel('Time (hrs)')
plt.ylabel('Fraction of Cells by Type')
plt.legend(loc=2)
plt.grid()
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.savefig('Fig.Cells.pdf',transparent=True)
plt.clf()

Ve = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    # plt.plot(f['Time'], f['Ve'], color='#666699',linewidth = 4.0, alpha=0.1)
    Ve[:, r - 1] = f['Ve']
plt.plot(f['Time'],np.mean(Ve,1),color='#666699',linewidth = 3.0,label='Viral Load')
plt.fill_between(f['Time'], np.amin(Ve, 1), np.amax(Ve, 1), facecolor='#666699', interpolate=True, alpha=0.2)
plt.title('Viral Load')
plt.xlabel('Time (hrs)')
plt.ylabel(r'log$_{10}$ PFU/ml')
plt.yscale('log')
# plt.ylim([1,None])
plt.legend(loc=2)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.ViralLoad.pdf',transparent=True)
plt.clf()

IFNe = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
    # plt.plot(f['Time'], f['IFNe'], color='#993300',linewidth = 4.0, alpha=0.1)
    IFNe[:, r - 1] = f['IFNe']
plt.plot(f['Time'],np.mean(IFNe,1),color='#993300',linewidth = 3.0,label='IFNe')
plt.fill_between(f['Time'], np.amin(IFNe, 1), np.amax(IFNe, 1), facecolor='#993300', interpolate=True, alpha=0.2)
plt.title('Extracellular IFN')
plt.xlabel('Time (hrs)')
plt.ylabel(r'[IFNe] $\mu$M')
plt.yscale('log')
plt.legend(loc=2)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.grid()
plt.savefig('Fig.TotalIFNe.pdf',transparent=True)
