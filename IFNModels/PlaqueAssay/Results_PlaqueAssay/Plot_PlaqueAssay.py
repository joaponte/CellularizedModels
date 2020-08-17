import numpy as np
import matplotlib.pyplot as plt


ts = 300
replicates = 20
f = np.genfromtxt('PlaqueAssay_%i.txt' % (1), skip_header=1, delimiter=',',max_rows = ts)
avgI1rd = np.zeros((ts, replicates))
avgI2rd = np.zeros((ts, replicates))
avgDrd = np.zeros((ts, replicates))
Beff = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('PlaqueAssay_%i.txt' % (r), skip_header=1, delimiter=',', max_rows = ts)
    plt.plot(f[:,0], f[:,1], color='orange',linewidth = 3.0, alpha=0.075)
    plt.plot(f[:,0], f[:,2], color='red',linewidth = 3.0, alpha=0.075)
    plt.plot(f[:,0], f[:,3], color='purple',linewidth = 3.0, alpha=0.075)
    avgI1rd[:, r - 1] = f[:,1]
    avgI2rd[:, r - 1] = f[:,2]
    avgDrd[:, r - 1] = f[:,3]
plt.plot(f[:,0],np.mean(avgI1rd,1),color='orange',linewidth = 3.0,label='I1')
plt.plot(f[:,0],np.mean(avgI2rd,1),color='red',linewidth = 3.0,label='I2')
plt.plot(f[:,0],np.mean(avgDrd,1),color='purple',linewidth = 3.0,label='D')
plt.title('Plaque Radius')
plt.savefig('Fig.Radius.pdf')
# plt.legend(loc=1)
# plt.savefig('Fig.Cells.pdf')
# plt.clf()
plt.show()
# Ve = np.zeros((ts, replicates))
# for r in range(1,replicates+1):
#     f = np.genfromtxt('FullModelCellular_%i.txt' % (r), skip_header=1, delimiter=',', names=CellNames,max_rows = ts)
#     plt.plot(f['Time'], f['Ve'], color='black',linewidth = 4.0, alpha=0.1)
#     Ve[:, r - 1] = f['Ve']
# plt.plot(f['Time'],np.mean(Ve,1),color='black',linewidth = 3.0,label='Viral Load')
# plt.legend(loc=2)
# plt.savefig('Fig.ViralLoad.pdf')
