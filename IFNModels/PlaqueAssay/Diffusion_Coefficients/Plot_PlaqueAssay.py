import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
ts = 468
replicates = 20
parameter = 'dcs'
p_I = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
p_V = [1.00, 2.00, 3.00, 4.00]
cell_diameter = 3.0
for i in range(len(p_I)):
    for j in range(len(p_V)):
        avgI1rd = np.zeros((ts, replicates))
        avgI2rd = np.zeros((ts, replicates))
        avgDrd = np.zeros((ts, replicates))
        # plt.title('Plaque Growth\n%s *= %.2f,\n%s *= %.2f' % ('IFN dc', p_I[i],'Virus dc', p_V[j]))
        for r in range(1,replicates+1):
            f = np.genfromtxt('PlaqueAssay_%s_%.2f_%.2f_%i.txt' % (parameter,p_I[i],p_V[j],r), skip_header=1, delimiter=',', names=names,max_rows=ts)
            avgI1rd[:, r - 1] = f['avgI1rd']
            avgI2rd[:, r - 1] = f['avgI2rd']
            avgDrd[:, r - 1] = f['avgDrd']
        plt.plot(f['Time'], np.mean(avgI1rd, 1) / cell_diameter, color='orange', linewidth=4.0, label='I1')
        plt.fill_between(f['Time'], np.amin(avgI1rd, 1) / cell_diameter, np.amax(avgI1rd, 1) / cell_diameter,
                         facecolor='orange', interpolate=True, alpha=0.2)
        plt.plot(f['Time'], np.mean(avgI2rd, 1) / cell_diameter, color='red', linewidth=4.0, label='I2')
        plt.fill_between(f['Time'], np.amin(avgI2rd, 1) / cell_diameter, np.amax(avgI2rd, 1) / cell_diameter,
                         facecolor='red', interpolate=True, alpha=0.2)
        plt.plot(f['Time'], np.mean(avgDrd, 1) / cell_diameter, color='purple', linewidth=4.0, label='D')
        plt.fill_between(f['Time'], np.amin(avgDrd, 1) / cell_diameter, np.amax(avgDrd, 1) / cell_diameter,
                         facecolor='purple', interpolate=True, alpha=0.2)
        plt.grid()
        plt.ylim([0, 100])
        plt.xlim([0, 80])
        # plt.legend(loc=2)
        if j == 0:
            plt.ylabel('Cell Diameters')
        if j > 0:
            plt.gca().set_yticklabels([])
        if i == 0:
            plt.xlabel('Time (hrs)')
        if i > 0:
            plt.gca().set_xticklabels([])
        plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
        plt.savefig('%s_%.2f_%.2f.pdf' % (parameter, p_I[i],p_V[j]), transparent = True)
        plt.clf()
