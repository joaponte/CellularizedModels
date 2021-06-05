import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 15})

# Baseline Slope
names = ['time', 'avgI1rd', 'avgI2rd', 'avgDrd', 'Beta', 'Beff']
t = 480
replicates = 20
DATA = np.zeros((19,4))
for j in range(2,replicates+1):
    MI1 = np.zeros((t, j))
    msI1 = np.zeros(j)
    for r in range(1, j+1):
        f = np.genfromtxt('Data/PlaqueAssay_%s_%.2f_%i.txt' % ('beta', 1.0, r), skip_header=1, delimiter=',',
                              names=names, max_rows=t)
        MI1[:, r - 1] = f['avgI1rd']
        # Calculate Slope AT THE END
        [mI1, y0I1] = np.polyfit(f['time'][-10:], f['avgI1rd'][-10:], 1)
        msI1[r - 1] = mI1
    Slope = np.mean(msI1)
    Stdev = np.std(msI1)
    Sterr = stats.sem(msI1)
    DATA[ j - 2,0] = j
    DATA[j - 2, 1] = Slope
    DATA[j - 2, 2] = Stdev
    DATA[j - 2, 3] = Sterr
plt.plot(DATA[:,0],DATA[:,2],label='SD',linewidth = 3.0)
plt.plot(DATA[:,0],DATA[:,3],label = 'SE',linewidth = 3.0)
plt.title('SD and SE of Plaque Growth Rate')
plt.ylabel('Magnitude')
plt.xlabel('Number of replicates')
plt.grid()
plt.gca().set_xticks([4,8,12,16,20])
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.savefig('Fig.Growth_Rate_SDSE.png', transparent = True)
plt.clf()

# Baseline Area Under the Curves
names = ['time', 'U', 'I1', 'I2', 'D', 'Ve', 'IFNe']
DATA_VeAUC = np.zeros((19,4))
DATA_VeMAX = np.zeros((19,4))
DATA_IFNeAUC = np.zeros((19,4))
DATA_IFNeMAX = np.zeros((19,4))
for j in range(2,replicates+1):
    rVeAUC = np.zeros(j)
    rIFNeAUC = np.zeros(j)
    rVeMAX = np.zeros(j)
    rIFNeMAX = np.zeros(j)
    for r in range(1, j+1):
        f = np.genfromtxt('Data/FullModelCellular_%s_%.2f_%i.txt' % ('beta', 1.0, r), skip_header=1, delimiter=',',
                          names=names, max_rows=t)
        rVeAUC[r - 1] = np.sum(f['Ve'])
        rVeMAX[r - 1] = np.amax(f['Ve'])
        rIFNeAUC[r - 1] = np.sum(f['IFNe'])
        rIFNeMAX[r - 1] = np.amax(f['IFNe'])
    VeAUC = np.mean(rVeAUC)
    VeAUCstd = np.std(rVeAUC)
    VeAUCste = stats.sem(rVeAUC)
    DATA_VeAUC[ j - 2,0] = j
    DATA_VeAUC[j - 2, 1] = VeAUC
    DATA_VeAUC[j - 2, 2] = VeAUCstd
    DATA_VeAUC[j - 2, 3] = VeAUCste
    VeMAX = np.mean(rVeMAX)
    VeMAXstd = np.std(rVeMAX)
    VeMAXste = stats.sem(rVeMAX)
    DATA_VeMAX[ j - 2,0] = j
    DATA_VeMAX[j - 2, 1] = VeMAX
    DATA_VeMAX[j - 2, 2] = VeMAXstd
    DATA_VeMAX[j - 2, 3] = VeMAXste
    IFNeAUC = np.mean(rIFNeAUC)
    IFNeAUCstd = np.std(rIFNeAUC)
    IFNeAUCste = stats.sem(rIFNeAUC)
    DATA_IFNeAUC[ j - 2,0] = j
    DATA_IFNeAUC[j - 2, 1] = IFNeAUC
    DATA_IFNeAUC[j - 2, 2] = IFNeAUCstd
    DATA_IFNeAUC[j - 2, 3] = IFNeAUCste
    IFNeMAX = np.mean(rIFNeMAX)
    IFNeMAXstd = np.std(rIFNeMAX)
    IFNeMAXste = stats.sem(rIFNeMAX)
    DATA_IFNeMAX[ j - 2,0] = j
    DATA_IFNeMAX[j - 2, 1] = IFNeMAX
    DATA_IFNeMAX[j - 2, 2] = IFNeMAXstd
    DATA_IFNeMAX[j - 2, 3] = IFNeMAXste

plt.plot(DATA_VeAUC[:,0],DATA_VeAUC[:,2],label='SD',linewidth = 3.0)
plt.plot(DATA_VeAUC[:,0],DATA_VeAUC[:,3],label = 'SE',linewidth = 3.0)
plt.title('SD and SE of VeAUC')
plt.ylabel('Magnitude')
plt.xlabel('Number of replicates')
plt.grid()
plt.gca().set_xticks([4,8,12,16,20])
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.savefig('Fig.VeAUC_SDSE.png', transparent = True)
plt.clf()

plt.plot(DATA_VeMAX[:,0],DATA_VeMAX[:,2],label='SD',linewidth = 3.0)
plt.plot(DATA_VeMAX[:,0],DATA_VeMAX[:,3],label = 'SE',linewidth = 3.0)
plt.title('SD and SE of VeMAX')
plt.ylabel('Magnitude')
plt.xlabel('Number of replicates')
plt.grid()
plt.gca().set_xticks([4,8,12,16,20])
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.savefig('Fig.VeMAX_SDSE.png', transparent = True)
plt.clf()

plt.plot(DATA_IFNeAUC[:,0],DATA_IFNeAUC[:,2],label='SD',linewidth = 3.0)
plt.plot(DATA_IFNeAUC[:,0],DATA_IFNeAUC[:,3],label = 'SE',linewidth = 3.0)
plt.title('SD and SE of IFNeAUC')
plt.ylabel('Magnitude')
plt.xlabel('Number of replicates')
plt.grid()
plt.gca().set_xticks([4,8,12,16,20])
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)

plt.savefig('Fig.IFNeAUC_SDSE.png', transparent = True)
plt.clf()
plt.plot(DATA_IFNeMAX[:,0],DATA_IFNeMAX[:,2],label='SD',linewidth = 3.0)
plt.plot(DATA_IFNeMAX[:,0],DATA_IFNeMAX[:,3],label = 'SE',linewidth = 3.0)
plt.title('SD and SE of IFNeMAX')
plt.ylabel('Magnitude')
plt.xlabel('Number of replicates')
plt.grid()
plt.gca().set_xticks([4,8,12,16,20])
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.savefig('Fig.IFNeMAX_SDSE.png', transparent = True)
plt.clf()