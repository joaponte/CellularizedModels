import numpy as np
import matplotlib.pyplot as plt

names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
t = 480
replicate = 20
parameter = 'k31'
p = [0.00, 1.00]
bpI1 = np.zeros(len(p))
bpI2 = np.zeros(len(p))
bpD= np.zeros(len(p))
for i in range(len(p)):
    MI1 = np.zeros((t, replicate))
    MI2 = np.zeros((t, replicate))
    MD = np.zeros((t, replicate))
    plt.title('Plaque Growth\n%s *= %.2f' % (parameter, p[i]))
    for r in range(1,replicate+1):
        f = np.genfromtxt('PlaqueAssay_%s_%.2f_%i.txt' % (parameter,p[i],r), skip_header=1, delimiter=',', names=names,max_rows=t)
        plt.plot(f['Time'],f['avgI1rd'],color='orange',linewidth = 3.0,alpha=0.075)
        plt.plot(f['Time'],f['avgI2rd'],color='red',linewidth = 3.0,alpha=0.075)
        plt.plot(f['Time'],f['avgDrd'],color='purple',linewidth = 3.0, alpha=0.075)
        MI1[:,r-1] = f['avgI1rd']
        MI2[:,r-1] = f['avgI2rd']
        MD[:,r-1] = f['avgDrd']
    plt.plot(f['Time'],np.mean(MI1,1),color='orange',linewidth = 3.0,label='I1')
    plt.plot(f['Time'],np.mean(MI2,1),color='red',linewidth = 3.0,label='I2')
    plt.plot(f['Time'],np.mean(MD,1),color='purple',linewidth = 3.0,label='D')
    plt.ylabel('Cell Diameters')
    plt.xlabel('Time (hrs)')
    plt.ylim([0, 150])
    plt.xlim([0, 80])
    plt.legend(loc=2)
    plt.savefig('%s_%.2f.pdf' % (parameter, p[i]))
    plt.clf()
    # Calculate Slope AT THE END
    [mI1,y0I1] = np.polyfit(f['Time'][-10:],np.mean(MI1,1)[-10:],1)
    bpI1[i] = mI1


plt.plot(p,bpI1,color='black',linewidth = 3.0, label = 'Plaque Growth Rate')
plt.title('Plaque Growth Rate')
plt.ylabel('Plaque Growth Rate (cell diameter/hrs)')
plt.xlabel('Parameter Multiplier')
plt.legend(loc=2)
plt.savefig('Growth_Rate_%s.pdf' % parameter)
