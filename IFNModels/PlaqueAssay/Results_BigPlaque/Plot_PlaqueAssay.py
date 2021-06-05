import numpy as np
import matplotlib.pyplot as plt

ts = 1000
replicates = 20
cell_diameter = 3.0
Names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
avgI1rd = np.zeros((ts, replicates))
avgI2rd = np.zeros((ts, replicates))
avgDrd = np.zeros((ts, replicates))
Beff = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    print('PlaqueAssay_%i.txt' % (r))
    f = np.genfromtxt('PlaqueAssay_%i.txt' % (r), skip_header=1, delimiter=',', names = Names, max_rows = ts)
    plt.plot(f['Time'], f['avgI1rd']/cell_diameter, color='orange',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['avgI2rd']/cell_diameter, color='red',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['avgDrd']/cell_diameter, color='purple',linewidth = 3.0, alpha=0.075)
    avgI1rd[:, r - 1] = f['avgI1rd']
    avgI2rd[:, r - 1] = f['avgI2rd']
    avgDrd[:, r - 1] = f['avgDrd']
plt.plot(f['Time'],np.mean(avgI1rd,1)/cell_diameter,color='orange',linewidth = 3.0,label='I1')
plt.plot(f['Time'],np.mean(avgI2rd,1)/cell_diameter,color='red',linewidth = 3.0,label='I2')
plt.plot(f['Time'],np.mean(avgDrd,1)/cell_diameter,color='purple',linewidth = 3.0,label='D')
plt.title('Plaque Radius')
plt.xlabel('Time (hrs)')
plt.ylabel('Cell Diameters')
plt.legend(loc=2)
plt.savefig('Fig.Radius.pdf')
plt.clf()

