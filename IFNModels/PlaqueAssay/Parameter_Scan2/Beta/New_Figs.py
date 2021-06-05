import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
plt.rcParams.update({'font.size': 15})

parameter = 'k31'

p = np.array([1.0,5.0,10.0,50.0,100.0])
VeAUC = np.array([4.96E+09,4.68E+09,6.19E+08,1.30E+07,3.21E+06])
VeAUCMin = np.array([7.93E+08,2.23E+08,5.00E+07,1.14E+05,2.77E+05])
VeAUCMax = np.array([5.34E+09,5.26E+09,8.10E+08,2.43E+07,5.78E+06])

plt.plot(p,VeAUC,color='#666699',linewidth = 4.0, label = 'Virus AUC')
plt.fill_between(p, VeAUCMin,VeAUCMax, facecolor='#666699', interpolate=True, alpha=0.2)
plt.title('Virus AUC')
plt.ylabel(r'Virus AUC (PFU*hr/mL)')
plt.xlabel('Parameter Multiplier')
plt.xscale('log')
# plt.yscale('log')
plt.grid()
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.gca().set_xticks(p)
plt.legend(loc=2)
plt.savefig('Fig.Virus_AUC_%s.pdf' % parameter,transparent=True)
plt.clf()

IFNeAUC = np.array([8.63E+05,5.39E+07,1.38E+09,8.08E+10,1.52E+11])
IFNeAUCMin = np.array([4.70E+02,5.32E+04,1.19E+09,8.05E+10,1.52E+11])
IFNeAUCMax = np.array([1.73E+07,9.46E+08,2.16E+09,8.17E+10,1.52E+11])

plt.plot(p,IFNeAUC,color='#993300',linewidth = 4.0, label = 'IFN AUC')
plt.fill_between(p, IFNeAUCMin,IFNeAUCMax, facecolor='#993300', interpolate=True, alpha=0.2)
plt.title('IFN AUC')
plt.ylabel(r'IFN AUC ($\mu$M*hr)')
plt.xlabel('Parameter Multiplier')
plt.xscale('log')
# plt.yscale('log')
plt.grid()
plt.gca().set_aspect(1.0 / plt.gca().get_data_ratio() * 1.0)
plt.gca().set_xticks(p)
plt.legend(loc=2)
plt.savefig('Fig.IFNe_AUC_%s.pdf' % parameter,transparent=True)
plt.clf()