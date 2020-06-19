import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

names2 = ['Time','T','I1A','I2A','DA','VA','I1B', 'I2B', 'DB', 'VB']
f2 = np.genfromtxt('V2CoinfectionODE_output.txt', skip_header=1, delimiter=',', names=names2)
names2 = ['Time','U','I1A','I2A','DA','VA','I1B', 'I2B', 'DB', 'VB']
f3 = np.genfromtxt('V2CoinfectionCellular_output.txt', skip_header=1, delimiter=',', names=names2)

T0 = f2['T'][0]
U0 = f3['U'][0]

days = 4.0
if f2['Time'][-1] < days:
    days = f2['Time'][-1]

plt.figure(figsize=(9.0, 4.5))
plt.subplot(1, 2, 1)
plt.plot(f2['Time'], f2['T'] / T0, '--', linewidth=2.0, color='blue')
plt.plot(f2['Time'], f2['I1A'] / T0, '--', linewidth=2.0, color='orange')
plt.plot(f2['Time'], f2['I2A'] / T0, '--', linewidth=2.0, color='red')
plt.plot(f2['Time'], f2['DA'] / T0, '--', linewidth=2.0, color='purple')
# plt.plot(f2['Time'], f2['I1B'] / T0, '--', linewidth=2.0, color='orange')
# plt.plot(f2['Time'], f2['I2B'] / T0, '--', linewidth=2.0, color='red')
# plt.plot(f2['Time'], f2['DB'] / T0, '--', linewidth=2.0, color='purple')
plt.plot(f3['Time'], f3['U'] / U0, linewidth=2.0, color="blue")
plt.plot(f3['Time'], f3['I1A'] / U0, linewidth=2.0, color='orange')
plt.plot(f3['Time'], f3['I2A'] / U0, linewidth=2.0, color='red')
plt.plot(f3['Time'], f3['DA'] / U0, linewidth=2.0, color='purple')
plt.plot(f3['Time'], f3['I1B'] / U0, linewidth=2.0, color='green')
plt.plot(f3['Time'], f3['I2B'] / U0, linewidth=2.0, color='magenta')
plt.plot(f3['Time'], f3['DB'] / U0, linewidth=2.0, color='brown')
plt.ylabel('Fraction of Cells')
plt.xlabel('Days')
plt.xlim([0, days])

l1 = Line2D([0], [0], linestyle='dashed', color='black', linewidth=2.0)
l2 = Line2D([0], [0], color='black', linewidth=2.0)
custom_lines = [l1, l2]
plt.legend(custom_lines, ['ODE model', 'Cellular model'], loc=1)

plt.subplot(1, 2, 2)
plt.plot(f2['Time'], f2['VB'], '--', linewidth=2.0, color='black')
plt.plot(f2['Time'], f2['VA'], '--', linewidth=2.0, color='black')
plt.plot(f3['Time'], f3['VA'], linewidth=2.0, color='black')
plt.plot(f3['Time'], f3['VB'], linewidth=2.0, color='black')

plt.yscale('log')
plt.xlim([0, days])
plt.ylim([10 ** 1, 10 ** 7])
plt.legend(['ODE Model', 'Cellular Model'], loc=1)
plt.ylabel('log TCID')
plt.xlabel('Days')

plt.tight_layout()
plt.savefig('V2TimeSeries.pdf')
plt.show()