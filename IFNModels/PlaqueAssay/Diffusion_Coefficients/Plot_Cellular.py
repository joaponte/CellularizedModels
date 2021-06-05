import numpy as np
import matplotlib.pyplot as plt

names = ['Time','U','I1','I2','D','Ve','IFNe']
t = 480
replicate = 20
parameter = 'IFNe_dc'
p = [1.00, 5.00, 10.00, 15.00]
Final_I1 = np.zeros(len(p))
Final_I2 = np.zeros(len(p))
Final_D = np.zeros(len(p))
Total_Ve = np.zeros(len(p))
Total_IFNe = np.zeros(len(p))
for i in range(len(p)):
    MI1 = np.zeros((1, replicate))
    MI2 = np.zeros((1, replicate))
    MD = np.zeros((1, replicate))
    TVe = np.zeros((1, replicate))
    TIFNe = np.zeros((1, replicate))
    for r in range(1,replicate+1):
        f = np.genfromtxt('FullModelCellular_%s_%.2f_%i.txt' % (parameter,p[i],r), skip_header=1, delimiter=',', names=names)
        MI1[:,r-1] = f['I1'][-1]
        MI2[:,r-1] = f['I2'][-1]
        MD[:,r-1] = f['D'][-1]
        TVe[:,r-1] = np.sum(f['Ve'])
        TIFNe[:, r - 1] = np.sum(f['IFNe'])
    Final_I1[i] = np.mean(MI1,1)
    Final_I2[i] = np.mean(MI2, 1)
    Final_D[i] = np.mean(MD, 1)
    Total_Ve[i] = np.mean(TVe, 1)
    Total_IFNe[i] = np.mean(TIFNe, 1)

plt.plot(p,Final_I1,color='orange',linewidth = 3.0, label = 'I1')
plt.plot(p,Final_I2,color='red',linewidth = 3.0, label = 'I2')
plt.plot(p,Final_D,color='purple',linewidth = 3.0,label = 'D')
plt.title('Final Fraction of Cells')
plt.ylabel('Final Fraction of Cells')
plt.xlabel('Parameter Multiplier')
plt.legend(loc=2)
plt.savefig('Final_Cells_%s.pdf' % parameter)
plt.clf()

plt.plot(p,Total_Ve,color='#666699',linewidth = 3.0, label = 'Extracellular Virus')
plt.title('Total Extracellular Virus')
plt.ylabel(r'Total Extracellular Virus (PFU/ml)')
plt.xlabel('Parameter Multiplier')
plt.legend(loc=2)
plt.savefig('Total_Virus_%s.pdf' % parameter)
plt.clf()

plt.plot(p,Total_IFNe,color='#993300',linewidth = 3.0, label = 'Extracellular IFN')
plt.title('Total Extracellular IFN')
plt.ylabel(r'Total Extracellular IFN ($\mu$M)')
plt.xlabel('Parameter Multiplier')
plt.legend(loc=2)
plt.savefig('Total_IFNe_%s.pdf' % parameter)