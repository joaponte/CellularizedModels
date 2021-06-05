import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parameters = ['beta', 'c', 'k', 'k12', 'k13', 'k14', 'k21', 'k31', 'k32', 'k33', 'k41', 'k42', 'k51', 'k61', 'k71',
              'k72', 'k73', 'n', 't3', 't4', 't5', 'vdc', 'idc']

# Baseline Slope
names = ['time', 'avgI1rd', 'avgI2rd', 'avgDrd', 'Beta', 'Beff']
t = 480
replicates = 20
MI1 = np.zeros((t, replicates))
MI2 = np.zeros((t, replicates))
MD = np.zeros((t, replicates))
for r in range(1, replicates + 1):
    f = np.genfromtxt('Data/PlaqueAssay_%s_%.2f_%i.txt' % ('beta', 1.0, r), skip_header=1, delimiter=',',
                      names=names, max_rows=t)
    MI1[:, r - 1] = f['avgI1rd']
    MI2[:, r - 1] = f['avgI2rd']
    MD[:, r - 1] = f['avgDrd']
# Calculate Slope AT THE END
[mI1, y0I1] = np.polyfit(f['time'][-10:], np.mean(MI1, 1)[-10:], 1)
BaselineSlope = mI1
print(mI1)

# Baseline Area Under the Curves
names = ['time', 'U', 'I1', 'I2', 'D', 'Ve', 'IFNe']
Ve = np.zeros((t, replicates))
IFNe = np.zeros((t, replicates))
for r in range(1, replicates + 1):
    f = np.genfromtxt('Data/FullModelCellular_%s_%.2f_%i.txt' % ('beta', 1.0, r), skip_header=1, delimiter=',',
                      names=names, max_rows=t)
    Ve[:, r - 1] = f['Ve']
    IFNe[:, r - 1] = f['IFNe']
BaselineVeAUC = np.sum(np.mean(Ve, 1))
BaselineIFNeAUC = np.sum(np.mean(IFNe, 1))
print(BaselineVeAUC)
print(BaselineIFNeAUC)

# Slopes
Outputs = ['Growth Rate']
names = ['Time', 'avgI1rd', 'avgI2rd', 'avgDrd', 'Beta', 'Beff']
Slope = np.zeros((len(parameters), 1))
for j in range(len(parameters)):
    p = [0.75, 1.25]
    bpI1 = np.zeros(len(p))
    for i in range(len(p)):
        MI1 = np.zeros((t, replicates))
        MI2 = np.zeros((t, replicates))
        MD = np.zeros((t, replicates))
        for r in range(1, replicates + 1):
            f = np.genfromtxt('Data/PlaqueAssay_%s_%.2f_%i.txt' % (parameters[j], p[i], r), skip_header=1,
                              delimiter=',', names=names, max_rows=t)
            MI1[:, r - 1] = f['avgI1rd']
            MI2[:, r - 1] = f['avgI2rd']
            MD[:, r - 1] = f['avgDrd']
        # Calculate Slope AT THE END
        [mI1, y0I1] = np.polyfit(f['Time'][-10:], np.mean(MI1, 1)[-10:], 1)
        bpI1[i] = mI1
    Slope[j] = max(abs(bpI1 - BaselineSlope)) / BaselineSlope

sns.heatmap(Slope, yticklabels=parameters, xticklabels=Outputs,square=True)
plt.savefig('Fig.Growth_Rate.png')
plt.show()
plt.clf()

# Area Under the Curves
# Outputs = ['Virus AUC','IFN AUC']
names = ['time', 'U', 'I1', 'I2', 'D', 'Ve', 'IFNe']
VeAUC = np.zeros((len(parameters), 1))
IFNeAUC = np.zeros((len(parameters), 1))
for j in range(len(parameters)):
    p = [0.75, 1.25]
    pVeAUC = np.zeros(len(p))
    pIFNeAUC = np.zeros(len(p))
    for i in range(len(p)):
        Ve = np.zeros((t, replicates))
        IFNe = np.zeros((t, replicates))
        for r in range(1, replicates + 1):
            f = np.genfromtxt('Data/FullModelCellular_%s_%.2f_%i.txt' % (parameters[j], p[i], r), skip_header=1,
                              delimiter=',', names=names, max_rows=t)
            Ve[:, r - 1] = f['Ve']
            IFNe[:, r - 1] = f['IFNe']
        pVeAUC[i] = np.sum(np.mean(Ve, 1))
        pIFNeAUC[i] = np.sum(np.mean(IFNe, 1))
    VeAUC[j] = max(abs(pVeAUC - BaselineVeAUC)) / BaselineVeAUC
    IFNeAUC[j] = max(abs(pIFNeAUC - BaselineIFNeAUC)) / BaselineIFNeAUC
Outputs = ['Virus AUC']
sns.heatmap(VeAUC, yticklabels=parameters, xticklabels=Outputs,square=True)
plt.savefig('Fig.VeAUC.png')
plt.show()
plt.clf()
Outputs = ['IFN AUC']
sns.heatmap(IFNeAUC, yticklabels=parameters, xticklabels=Outputs,square=True)
plt.savefig('Fig.IFNeAUC.png')
plt.show()
plt.clf()

Outputs = ['Growth Rate','Virus AUC','IFN AUC']
heatmap_matrix = np.zeros((3,len(parameters)))
heatmap_matrix[0,:,] = np.log10(Slope[:,0])
heatmap_matrix[1,:,] = np.log10(VeAUC[:,0])
heatmap_matrix[2,:,] = np.log10(IFNeAUC[:,0])
# heatmap_matrix[0,:,] = Slope[:,0]
# heatmap_matrix[1,:,] = VeAUC[:,0]
# heatmap_matrix[2,:,] = IFNeAUC[:,0]
sns.heatmap(heatmap_matrix.transpose(), yticklabels=parameters, xticklabels=Outputs,square=True)
plt.savefig('Fig.Panel.png')
plt.show()
plt.clf()

# cmap=ListedColormap(['white'])
# fmt='g'