import numpy as np
import matplotlib.pyplot as plt


ts = 300
replicates = 20
Names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
avgI1rd = np.zeros((ts, replicates))
avgI2rd = np.zeros((ts, replicates))
avgDrd = np.zeros((ts, replicates))
Beff = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('PlaqueAssay_%i.txt' % (r), skip_header=1, delimiter=',', names = Names, max_rows = ts)
    plt.plot(f['Time'], f['avgI1rd'], color='orange',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['avgI2rd'], color='red',linewidth = 3.0, alpha=0.075)
    plt.plot(f['Time'], f['avgDrd'], color='purple',linewidth = 3.0, alpha=0.075)
    avgI1rd[:, r - 1] = f['avgI1rd']
    avgI2rd[:, r - 1] = f['avgI2rd']
    avgDrd[:, r - 1] = f['avgDrd']
plt.plot(f['Time'],np.mean(avgI1rd,1),color='orange',linewidth = 3.0,label='I1')
plt.plot(f['Time'],np.mean(avgI2rd,1),color='red',linewidth = 3.0,label='I2')
plt.plot(f['Time'],np.mean(avgDrd,1),color='purple',linewidth = 3.0,label='D')
plt.title('Plaque Radius')

plt.legend(loc=1)
plt.savefig('Fig.Radius.pdf')
plt.clf()

Beff = np.zeros((ts, replicates))
for r in range(1,replicates+1):
    f = np.genfromtxt('PlaqueAssay_%i.txt' % (r), skip_header=1, delimiter=',', names = Names, max_rows = ts)
    #plt.plot(f['Time'], f['Beff'], color='black',linewidth = 4.0, alpha=0.1)
    Beff[:, r - 1] = f['Beff']
plt.plot(f['Time'],np.mean(Beff,1),color='black',linewidth = 3.0,label='Effective Infectivity')
plt.plot(f['Time'],f['Beta'],'--',color='black',linewidth = 3.0,label='Theoretical Infectivity')
plt.legend(loc=2)
plt.ylim([0,max(f['Beta'])*1.1])
plt.show()
#plt.savefig('Fig.Beff.pdf')
