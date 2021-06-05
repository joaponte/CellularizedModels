import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

parameters = ['beta', 'c', 'k', 'k12', 'k13', 'k14', 'k21', 'k31', 'k32', 'k33', 'k41', 'k42', 'k51', 'k61', 'k71',
              'k72', 'k73', 'n', 't3', 't4', 't5', 'vdc', 'idc']

# Baseline Slope
names = ['time', 'avgI1rd', 'avgI2rd', 'avgDrd', 'Beta', 'Beff']
t = 480
replicates = 20
BaselinemsI1 = np.zeros(replicates)
for r in range(replicates):
    f = np.genfromtxt('Data/PlaqueAssay_%s_%.2f_%i.txt' % ('beta', 1.0, r), skip_header=1, delimiter=',',
                      names=names, max_rows=t)
    # Calculate Slope AT THE END
    [mI1, y0I1] = np.polyfit(f['time'][-10:], f['avgI1rd'][-10:], 1)
    BaselinemsI1[r] = mI1
BaselineSlope = np.mean(BaselinemsI1)
BaselineSlopestd = np.std(BaselinemsI1)
print('BaselineGrowthRate',BaselineSlope)
print('BaselineGrowthRateSTD',BaselineSlopestd)

# Baseline Area Under the Curves
names = ['time', 'U', 'I1', 'I2', 'D', 'Ve', 'IFNe']
rVeAUC = np.zeros(replicates)
rIFNeAUC = np.zeros(replicates)
rVeMAX = np.zeros(replicates)
rIFNeMAX = np.zeros(replicates)
for r in range(replicates):
    f = np.genfromtxt('Data/FullModelCellular_%s_%.2f_%i.txt' % ('beta', 1.0, r), skip_header=1, delimiter=',',
                      names=names, max_rows=t)
    rVeAUC[r] = np.sum(f['Ve'])
    rVeMAX[r] = np.amax(f['Ve'])
    rIFNeAUC[r] = np.sum(f['IFNe'])
    rIFNeMAX[r] = np.amax(f['IFNe'])
BaselineVeAUC = np.mean(rVeAUC)
BaselineVeAUCstd = np.std(rVeAUC)
BaselineVeMAX = np.mean(rVeMAX)
BaselineVeMAXstd = np.std(rVeMAX)
BaselineIFNeAUC = np.mean(rIFNeAUC)
BaselineIFNeAUCstd = np.std(rIFNeAUC)
BaselineIFNeMAX = np.mean(rIFNeMAX)
BaselineIFNeMAXstd = np.std(rIFNeMAX)
print('BaselineVeAUC',BaselineVeAUC)
print('BaselineVeAUCstd',BaselineVeAUCstd)
print('BaselineVeMAX',BaselineVeMAX)
print('BaselineVeMAXstd',BaselineVeMAXstd)
print('BaselineIFNeAUC',BaselineIFNeAUC)
print('BaselineIFNeAUCstd',BaselineIFNeAUCstd)
print('BaselineIFNeMAX',BaselineIFNeMAX)
print('BaselineIFNeMAXstd',BaselineIFNeMAXstd)

# Slopes
names = ['time', 'avgI1rd', 'avgI2rd', 'avgDrd', 'Beta', 'Beff']
SlopeUp = np.zeros((len(parameters), 1))
SlopeUpPValue = np.zeros((len(parameters), 1))
SlopeDown = np.zeros((len(parameters), 1))
SlopeDownPValue = np.zeros((len(parameters), 1))
for j in range(len(parameters)):
    p = [0.75, 1.25]
    bpI1 = np.zeros(len(p))
    pval = np.zeros(len(p))
    for i in range(len(p)):
        msI1 = np.zeros(replicates)
        for r in range(replicates):
            f = np.genfromtxt('Data/PlaqueAssay_%s_%.2f_%i.txt' % (parameters[j], p[i], r), skip_header=1,
                              delimiter=',', names=names, max_rows=t)
            # Calculate Slope AT THE END
            [mI1, y0I1] = np.polyfit(f['time'][-10:], f['avgI1rd'][-10:], 1)
            msI1[r] = mI1
        bpI1[i] = np.mean(msI1)
        pval[i] = stats.ttest_ind(BaselinemsI1,msI1,equal_var=False).pvalue
    SlopeDown[j] = bpI1[0]
    SlopeDownPValue[j] = pval[0]
    SlopeUp[j] = bpI1[1]
    SlopeUpPValue[j] = pval[1]

Outputs = ['Baseline','stdev','Down','Diff_Down','p-value','Up','Diff_Up','p-value','Diff_Diff']
heatmap_matrix = np.zeros((9,len(parameters)))
heatmap_matrix[0,:,] = BaselineSlope
heatmap_matrix[1,:,] = BaselineSlopestd
heatmap_matrix[2,:,] = SlopeDown[:,0]
heatmap_matrix[3,:,] = BaselineSlope - SlopeDown[:,0]
heatmap_matrix[4,:,] = SlopeDownPValue[:,0]
heatmap_matrix[5,:,] = SlopeUp[:,0]
heatmap_matrix[6,:,] = SlopeUp[:,0] - BaselineSlope
heatmap_matrix[7,:,] = SlopeUpPValue[:,0]
heatmap_matrix[8,:,] = (SlopeUp[:,0] - 2* BaselineSlope + SlopeDown[:,0])/2.0
np.savetxt('Tab.k31.GrowthRate.txt',heatmap_matrix.transpose(),delimiter = ',')
sns.heatmap(heatmap_matrix.transpose(), yticklabels=parameters, xticklabels=Outputs, cmap = 'seismic',
            annot= True ,fmt='.3f',center=0.0)
plt.title('Growth Rate')
plt.tight_layout()
plt.savefig('Tab.Growth_Rate.png')
plt.clf()

# Area Under the Curves and Maximums
names = ['time', 'U', 'I1', 'I2', 'D', 'Ve', 'IFNe']
VeAUCup = np.zeros((len(parameters), 1))
VeAUCupPvalue = np.zeros((len(parameters), 1))
VeAUCdown = np.zeros((len(parameters), 1))
VeAUCdownPvalue = np.zeros((len(parameters), 1))
VeMAXup = np.zeros((len(parameters), 1))
VeMAXupPvalue = np.zeros((len(parameters), 1))
VeMAXdown = np.zeros((len(parameters), 1))
VeMAXdownPvalue = np.zeros((len(parameters), 1))
IFNeAUCup = np.zeros((len(parameters), 1))
IFNeAUCupPvalue = np.zeros((len(parameters), 1))
IFNeAUCdown = np.zeros((len(parameters), 1))
IFNeAUCdownPvalue = np.zeros((len(parameters), 1))
IFNeMAXup = np.zeros((len(parameters), 1))
IFNeMAXupPvalue = np.zeros((len(parameters), 1))
IFNeMAXdown = np.zeros((len(parameters), 1))
IFNeMAXdownPvalue = np.zeros((len(parameters), 1))
for j in range(len(parameters)):
    p = [0.75, 1.25]
    pVeAUC = np.zeros(len(p))
    pVeAUCPval = np.zeros(len(p))
    pVeMAX = np.zeros(len(p))
    pVeMAXPval = np.zeros(len(p))
    pIFNeAUC = np.zeros(len(p))
    pIFNeAUCPval = np.zeros(len(p))
    pIFNeMAX = np.zeros(len(p))
    pIFNeMAXPval = np.zeros(len(p))
    for i in range(len(p)):
        iVeAUC = np.zeros(replicates)
        iIFNeAUC = np.zeros(replicates)
        iVeMAX = np.zeros(replicates)
        iIFNeMAX = np.zeros(replicates)
        for r in range(replicates):
            f = np.genfromtxt('Data/FullModelCellular_%s_%.2f_%i.txt' % (parameters[j], p[i], r), skip_header=1,
                              delimiter=',', names=names, max_rows=t)
            iVeAUC[r] = np.sum(f['Ve'])
            iIFNeAUC[r] = np.sum(f['IFNe'])
            iVeMAX[r] = np.amax(f['Ve'])
            iIFNeMAX[r] = np.amax(f['IFNe'])
        pVeAUC[i] = np.mean(iVeAUC)
        pVeAUCPval[i] = stats.ttest_ind(rVeAUC,iVeAUC,equal_var=False).pvalue
        pVeMAX[i] = np.mean(iVeMAX)
        pVeMAXPval[i] = stats.ttest_ind(rVeMAX, iVeMAX,equal_var=False).pvalue
        pIFNeAUC[i] = np.mean(iIFNeAUC)
        pIFNeAUCPval[i] = stats.ttest_ind(rIFNeAUC, iIFNeAUC,equal_var=False).pvalue
        pIFNeMAX[i] = np.mean(iIFNeMAX)
        pIFNeMAXPval[i] = stats.ttest_ind(rIFNeMAX, iIFNeMAX,equal_var=False).pvalue
    VeAUCdown[j] = pVeAUC[0] / 1.0e10
    VeAUCdownPvalue[j] = pVeAUCPval[0]
    VeAUCup[j] = pVeAUC[1] / 1.0e10
    VeAUCupPvalue[j] = pVeAUCPval[1]
    VeMAXdown[j] = pVeMAX[1] / 1.0e8
    VeMAXdownPvalue[j] = pVeMAXPval[0]
    VeMAXup[j] = pVeMAX[0] / 1.0e8
    VeMAXupPvalue[j] = pVeMAXPval[1]
    IFNeAUCdown[j] = pIFNeAUC[0] / 1.0e10
    IFNeAUCdownPvalue[j] = pIFNeAUCPval[0]
    IFNeAUCup[j] = pIFNeAUC[1] / 1.0e10
    IFNeAUCupPvalue[j] = pIFNeAUCPval[1]
    IFNeMAXdown[j] = pIFNeMAX[0] / 1.0e8
    IFNeMAXdownPvalue[j] = pIFNeMAXPval[0]
    IFNeMAXup[j] = pIFNeMAX[1] / 1.0e8
    IFNeMAXupPvalue[j] = pIFNeMAXPval[1]

BaselineVeAUC /= 1.0e10
BaselineVeAUCstd /= 1.0e10
Outputs = ['Baseline','stdev','Down','Diff_Down','p-value','Up','Diff_Up','p-value','Diff_Diff']
heatmap_matrix = np.zeros((9,len(parameters)))
heatmap_matrix[0,:,] = BaselineVeAUC
heatmap_matrix[1,:,] = BaselineVeAUCstd
heatmap_matrix[2,:,] = VeAUCdown[:,0]
heatmap_matrix[3,:,] = BaselineVeAUC - VeAUCdown[:,0]
heatmap_matrix[4,:,] = VeAUCdownPvalue[:,0]
heatmap_matrix[5,:,] = VeAUCup[:,0]
heatmap_matrix[6,:,] = VeAUCup[:,0] - BaselineVeAUC
heatmap_matrix[7,:,] = VeAUCupPvalue[:,0]
heatmap_matrix[8,:,] = (VeAUCup[:,0] - 2 * BaselineVeAUC + VeAUCdown[:,0])/2.0
np.savetxt('Tab.k31.VeAUC.txt',heatmap_matrix.transpose(),delimiter = ',')
sns.heatmap(heatmap_matrix.transpose(), yticklabels=parameters, xticklabels=Outputs, cmap = 'seismic',
            annot= True ,fmt='.3f',center=0.0)
plt.title('Ve AUC')
plt.tight_layout()
plt.savefig('Tab.VeAUC.png')
plt.clf()

BaselineVeMAX /= 1.0e8
BaselineVeMAXstd /= 1.0e8
Outputs = ['Baseline','stdev','Down','Diff_Down','p-value','Up','Diff_Up','p-value','Diff_Diff']
heatmap_matrix = np.zeros((9,len(parameters)))
heatmap_matrix[0,:,] = BaselineVeMAX
heatmap_matrix[1,:,] = BaselineVeMAXstd
heatmap_matrix[2,:,] = VeMAXdown[:,0]
heatmap_matrix[3,:,] = BaselineVeMAX - VeMAXdown[:,0]
heatmap_matrix[4,:,] = VeMAXdownPvalue[:,0]
heatmap_matrix[5,:,] = VeMAXup[:,0]
heatmap_matrix[6,:,] = VeMAXup[:,0] - BaselineVeMAX
heatmap_matrix[7,:,] = VeMAXupPvalue[:,0]
heatmap_matrix[8,:,] = (VeMAXup[:,0] - 2 * BaselineVeMAX + VeMAXdown[:,0])/2.0
np.savetxt('Tab.k31.VeMAX.txt',heatmap_matrix.transpose(),delimiter = ',')
sns.heatmap(heatmap_matrix.transpose(), yticklabels=parameters, xticklabels=Outputs, cmap = 'seismic',
            annot= True ,fmt='.3f',center=0.0)
plt.title('Ve MAX')
plt.tight_layout()
plt.savefig('Tab.VeMAX.png')
plt.clf()

BaselineIFNeAUC /= 1.0e10
BaselineIFNeAUCstd /= 1.0e10
Outputs = ['Baseline','stdev','Down','Diff_Down','p-value','Up','Diff_Up','p-value','Diff_Diff']
heatmap_matrix = np.zeros((9,len(parameters)))
heatmap_matrix[0,:,] = BaselineIFNeAUC
heatmap_matrix[1,:,] = BaselineIFNeAUCstd
heatmap_matrix[2,:,] = IFNeAUCdown[:,0]
heatmap_matrix[3,:,] = BaselineIFNeAUC - IFNeAUCdown[:,0]
heatmap_matrix[4,:,] = IFNeAUCdownPvalue[:,0]
heatmap_matrix[5,:,] = IFNeAUCup[:,0]
heatmap_matrix[6,:,] = IFNeAUCup[:,0] - BaselineIFNeAUC
heatmap_matrix[7,:,] = IFNeAUCupPvalue[:,0]
heatmap_matrix[8,:,] = (IFNeAUCup[:,0] - 2 * BaselineIFNeAUC + IFNeAUCdown[:,0])/2.0
np.savetxt('Tab.k31.IFNeAUC.txt',heatmap_matrix.transpose(),delimiter = ',')
sns.heatmap(heatmap_matrix.transpose(), yticklabels=parameters, xticklabels=Outputs, cmap = 'seismic',
            annot= True ,fmt='.3f',center=0.0)
plt.title('IFNe AUC')
plt.tight_layout()
plt.savefig('Tab.IFNeAUC.png')
plt.clf()

BaselineIFNeMAX /= 1.0e8
BaselineIFNeMAXstd /= 1.0e8
Outputs = ['Baseline','stdev','Down','Diff_Down','p-value','Up','Diff_Up','p-value','Diff_Diff']
heatmap_matrix = np.zeros((9,len(parameters)))
heatmap_matrix[0,:,] = BaselineIFNeMAX
heatmap_matrix[1,:,] = BaselineIFNeMAXstd
heatmap_matrix[2,:,] = IFNeMAXdown[:,0]
heatmap_matrix[3,:,] = BaselineIFNeMAX - IFNeMAXdown[:,0]
heatmap_matrix[4,:,] = IFNeMAXdownPvalue[:,0]
heatmap_matrix[5,:,] = IFNeMAXup[:,0]
heatmap_matrix[6,:,] = IFNeMAXup[:,0] - BaselineIFNeMAX
heatmap_matrix[7,:,] = IFNeMAXupPvalue[:,0]
heatmap_matrix[8,:,] = (IFNeMAXup[:,0] - 2 * BaselineIFNeMAX + IFNeMAXdown[:,0])/2.0
np.savetxt('Tab.k31.IFNeMAX.txt',heatmap_matrix.transpose(),delimiter = ',')
sns.heatmap(heatmap_matrix.transpose(), yticklabels=parameters, xticklabels=Outputs, cmap = 'seismic',
            annot= True ,fmt='.3f',center=0.0)
plt.title('IFNe MAX')
plt.tight_layout()
plt.savefig('Tab.IFNeMAX.png')