import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

ts = 480
replicates = 20
cell_diameter = 3.0
Names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
avgI1rd = np.zeros((ts, replicates))
avgI2rd = np.zeros((ts, replicates))
avgDrd = np.zeros((ts, replicates))
Beff = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('PlaqueAssay_%i.txt' % (r), skip_header=1, delimiter=',', names = Names, max_rows = ts)
    # plt.plot(f['Time'], f['avgI1rd']/cell_diameter, color='orange',linewidth = 3.0, alpha=0.075)
    # plt.plot(f['Time'], f['avgI2rd']/cell_diameter, color='red',linewidth = 3.0, alpha=0.075)
    # plt.plot(f['Time'], f['avgDrd']/cell_diameter, color='purple',linewidth = 3.0, alpha=0.075)
    avgI1rd[:, r - 1] = f['avgI1rd']
    avgI2rd[:, r - 1] = f['avgI2rd']
    avgDrd[:, r - 1] = f['avgDrd']
plt.plot(f['Time'],np.mean(avgI1rd,1)/cell_diameter,color='orange',linewidth = 4.0,label='I1')
plt.fill_between(f['Time'], np.amin(avgI1rd, 1)/cell_diameter, np.amax(avgI1rd, 1)/cell_diameter, facecolor='orange', interpolate=True, alpha=0.2)
plt.plot(f['Time'],np.mean(avgI2rd,1)/cell_diameter,color='red',linewidth = 4.0,label='I2')
plt.fill_between(f['Time'], np.amin(avgI2rd, 1)/cell_diameter, np.amax(avgI2rd, 1)/cell_diameter, facecolor='red', interpolate=True, alpha=0.2)
plt.plot(f['Time'],np.mean(avgDrd,1)/cell_diameter,color='purple',linewidth = 4.0,label='D')
plt.fill_between(f['Time'], np.amin(avgDrd, 1)/cell_diameter, np.amax(avgDrd, 1)/cell_diameter, facecolor='purple', interpolate=True, alpha=0.2)
plt.grid()
plt.title('Plaque Radius')
plt.xlabel('Time (hrs)')
plt.ylabel('Cell Diameters')
plt.ylim([0,100])
plt.xlim([0,80])
plt.legend(loc=2)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.savefig('Fig.Radius.pdf', transparent = True)
plt.clf()

# Beff = np.zeros((ts, replicates))
# for r in range(1,replicates+1):
#     f = np.genfromtxt('PlaqueAssay_%i.txt' % (r), skip_header=1, delimiter=',', names = Names, max_rows = ts)
#     plt.plot(f['Time'][f['Beff']>0], f['Beff'][f['Beff']>0], color='black',linewidth = 3.0, alpha=0.075)
#     Beff[:, r - 1] = f['Beff']
# plt.plot(f['Time'][np.mean(Beff,1)>0],np.mean(Beff,1)[np.mean(Beff,1)>0],color='black',linewidth = 3.0,label='Effective Infectivity')
# plt.plot(f['Time'],f['Beta'],'--',color='black',linewidth = 3.0,label='Theoretical Infectivity')
# plt.legend(loc=1)
# plt.yscale('log')
# plt.savefig('Fig.Beff.pdf')

