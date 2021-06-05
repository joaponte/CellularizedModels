import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
plt.rcParams.update({'font.size': 15})

parameter = 'k11'

p = np.array([0.00,0.25,0.5,0.75,1.0])
VeAUC = np.array([4.89E+09,2.18E+09,1.33E+09,7.85E+08,5.65E+08])
VeAUCMin = np.array([2.71E+09,4.96E+08,8.67E+08,3.70E+08,3.47E+08])
VeAUCMax = np.array([5.36E+09,2.49E+09,1.45E+09,9.02E+08,6.11E+08])

plt.plot(p,VeAUC,color='#666699',linewidth = 4.0, label = 'Virus AUC')
plt.fill_between(p, VeAUCMin,VeAUCMax, facecolor='#666699', interpolate=True, alpha=0.2)
plt.title('Virus AUC')
plt.ylabel(r'Virus AUC (PFU*hr/mL)')
plt.xlabel('Parameter Multiplier')
plt.yscale('log')
plt.grid()
plt.gca().set_yticks([10**8,10**9,10**10])
plt.ylim([10**(7.9),10**(10.1)])
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.gca().set_xticks(p)
plt.legend(loc=2)
plt.savefig('Fig.Virus_AUC_%s.pdf' % parameter,transparent=True)
plt.show()

IFNeAUC = np.array([4.54E+03,3.82E+07,6.36E+07,8.17E+07,9.65E+07])
IFNeAUCMin = np.array([4.72E+02,3.47E+07,5.80E+07,7.28E+07,8.03E+07])
IFNeAUCMax = np.array([7.75E+04,6.60E+07,6.80E+07,8.78E+07,1.00E+08])

plt.plot(p,IFNeAUC,color='#993300',linewidth = 4.0, label = 'IFN AUC')
plt.fill_between(p, IFNeAUCMin,IFNeAUCMax, facecolor='#993300', interpolate=True, alpha=0.2)
plt.title('IFN AUC')
plt.ylabel(r'IFN AUC ($\mu$M*hr)')
plt.xlabel('Parameter Multiplier')
# plt.yscale('log')
plt.grid()
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.gca().set_xticks(p)
plt.legend(loc=2)
plt.savefig('Fig.IFNe_AUC_%s.pdf' % parameter,transparent=True)
plt.clf()