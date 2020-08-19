import numpy as np
import matplotlib.pyplot as plt

names = ['Time','avgI1rd','avgI2rd','avgDrd','Beta','Beff']
t = 300
replicate = 20
parameter = 'k71'
p = [0.1,0.5,1.0,2.0,10.0]
bpI1 = np.zeros(len(p))
bpI2 = np.zeros(len(p))
bpD= np.zeros(len(p))
for i in range(len(p)):
    MI1 = np.zeros((t, replicate))
    MI2 = np.zeros((t, replicate))
    MD = np.zeros((t, replicate))
    plt.subplot(2,3,i+1)
    plt.title(r'$%s = %.2f$' % (parameter, p[i]))
    for r in range(1,replicate+1):
        f = np.genfromtxt('PlaqueAssay_%s_%.2f_%i.txt' % (parameter,p[i],r), skip_header=1, delimiter=',', names=names,max_rows=t)
        plt.plot(f['Time'],f['avgI1rd'],color='orange',linewidth = 3.0,alpha=0.075)
        plt.plot(f['Time'],f['avgI2rd'],color='red',linewidth = 3.0,alpha=0.075)
        plt.plot(f['Time'],f['avgDrd'],color='purple',linewidth = 3.0, alpha=0.075)
        MI1[:,r-1] = f['avgI1rd']
        MI2[:,r-1] = f['avgI2rd']
        MD[:,r-1] = f['avgDrd']
    plt.plot(f['Time'],np.mean(MI1,1),color='orange',linewidth = 3.0)
    plt.plot(f['Time'],np.mean(MI2,1),color='red',linewidth = 3.0)
    plt.plot(f['Time'],np.mean(MD,1),color='purple',linewidth = 3.0)
    # mMI1 = np.mean(MI1,1)
    # [mI1,y0I1] = np.polyfit(f['Time'][np.nonzero(mMI1)],mMI1[np.nonzero(mMI1)],1)
    # bpI1[i] = mI1
    # mMI2 = np.mean(MI2,1)
    # [mI2,y0I2] = np.polyfit(f['Time'][np.nonzero(mMI2)],mMI2[np.nonzero(mMI2)],1)
    # bpI2[i] = mI2
    # mMD = np.mean(MD,1)
    # [mD,y0D] = np.polyfit(f['Time'][np.nonzero(mMD)],mMD[np.nonzero(mMD)],1)
    # bpD[i] = mD
    # x = np.arange(51)
    # yI1 = y0I1 + mI1 * x
    # yI2 = y0I2 + mI2 * x
    # yD = y0D + mD * x
    # plt.plot(x,yI1,color='orange',linewidth = 3.0)
    # plt.plot(x,yI2,color='red',linewidth = 3.0)
    # plt.plot(x,yD,color='purple',linewidth = 3.0)
    # plt.xlim([0,50])
    plt.ylim([0,50])
# plt.subplot(2,3,i+2)
# plt.plot(p,bpI1,color='orange',linewidth = 3.0)
# plt.plot(p,bpI2,color='red',linewidth = 3.0)
# plt.plot(p,bpD,color='purple',linewidth = 3.0)
# plt.xscale('log')
# plt.title("Plaque Growth")
plt.tight_layout()
plt.savefig('%s.pdf' % parameter)
plt.show()