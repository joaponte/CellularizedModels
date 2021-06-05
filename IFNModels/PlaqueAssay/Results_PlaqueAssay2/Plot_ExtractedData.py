import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

TimesB = [24, 36 ,48, 60, 72]
FigureB = [7.52688172, 27.15053763, 43.01075269, 64.24731183, 76.88172043]
TimesBD = [48, 60, 72]
FigureBD = [23.11827957,47.04301075,61.02150538]

[mFB,y0FB] = np.polyfit(TimesB,FigureB,1)
LineTimes = np.arange(0,80,0.1)
LineFigureB = LineTimes * mFB + y0FB

[mFBD,y0FBD] = np.polyfit(TimesBD,FigureBD,1)
LineFigureBD = LineTimes * mFBD + y0FBD

plt.plot(LineTimes,LineFigureB,'--',color='orange',linewidth = 4.0)
plt.plot(TimesB,FigureB,'s',color='black',markersize = 8.0,label='Infected')
plt.plot(LineTimes,LineFigureBD,'--' ,color='purple',linewidth = 4.0)
plt.plot(TimesBD,FigureBD,'o',color='black',markersize = 8.0,label='Dead')
plt.grid()
plt.title('Plaque Radius')
plt.xlabel('Time (hrs)')
plt.ylabel('Cell Diameters')
plt.ylim([0,100])
plt.xlim([0,80])
plt.legend(loc=2)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.savefig('FigB.Radius.pdf', transparent = True)
plt.clf()

TimesC = [36 ,48, 60, 72]
FigureC = [16.07142857, 27.38095238, 44.3452381, 48.80952381]
TimesCD = [48, 60, 72]
FigureCD = [10.11904762, 29.16666667, 36.9047619]

[mFC,y0FC] = np.polyfit(TimesC,FigureC,1)
LineTimes = np.arange(0,80,0.1)
LineFigureC = LineTimes * mFC + y0FC

[mFCD,y0FCD] = np.polyfit(TimesCD,FigureCD,1)
LineFigureCD = LineTimes * mFCD + y0FCD

plt.plot(LineTimes,LineFigureC,'--' ,color='orange',linewidth = 4.0)
plt.plot(TimesC,FigureC,'s',color='black',markersize = 8.0,label='Infected')
plt.plot(LineTimes,LineFigureCD,'--' ,color='purple',linewidth = 4.0)
plt.plot(TimesCD,FigureCD,'o',color='black',markersize = 8.0,label='Dead')
plt.grid()
plt.title('Plaque Radius')
plt.xlabel('Time (hrs)')
plt.ylabel('Cell Diameters')
plt.ylim([0,100])
plt.xlim([0,80])
plt.legend(loc=2)
plt.gca().set_aspect(1.0/plt.gca().get_data_ratio()*1.0)
plt.savefig('FigC.Radius.pdf', transparent = True)