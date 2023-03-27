import numpy as np
import matplotlib.pyplot as plt

BS = 1 # bootstrapping or not

proteins = ['PTEN_Pt','PTEN_C2','SHIP2_Pt','SHIP2_C2']
names = ['Ptase$_\mathrm{PTEN}$','C2$_\mathrm{PTEN}$','Ptase$_\mathrm{SHIP2}$','C2$_\mathrm{SHIP2}$']
colors = ['cyan','blue','magenta','red']
linewidth = 1.5
fontsize = 12.7
for (protein,name,color) in zip(proteins,names,colors):
    
    skip_first=10
    skip_last=4
    if BS:
        filename = '%s/umbrella/bsResult.xvg' % protein
        dist,E,dE = np.genfromtxt(filename,skip_header=18,unpack=True)
    else:
        filename = '%s/umbrella/profile.xvg' % protein
        dist,E = np.genfromtxt(filename,skip_header=17+skip_first,skip_footer=skip_last,unpack=True)
    zero = np.mean(E[-20:])
    E_zero = E-zero
    dPMF = np.amin(E_zero)
    #plt.plot(dist,E_zero,marker='*',color=color)
    #time_us = time/1000
    #plt.plot(dist,E_zero,color=color,label=r'%s mode %d: $\Delta$PMF = %1.0f kJ/mol,time = %d $\mu$s' % (protein_name,mode,dPMF,time_us))
    
    if BS:
        plt.fill_between(dist,E_zero-dE,E_zero+dE,color=color,alpha=0.5)
        #plt.plot(dist[-20:],E_zero[-20:],marker='.',color='green')
        idx_min=np.argmin(E_zero)
        #plt.plot(dist[idx],E_zero[idx],marker='.',color='green')
        dd_zero = np.sqrt(np.sum(dE[-20:]**2))/np.sqrt(20)
        dd_min  = dE[idx_min]
        ddPMF   = np.sqrt(dd_zero**2+dd_min**2)
        #plt.plot(dist,E_zero,linewidth=linewidth,color=color,label=r'%s, $\Delta$PMF = %1.0f +/- %1.0f kJ/mol' % (protein_name,dPMF,ddPMF))
        plt.plot(dist,E_zero,linewidth=linewidth,color=color,label=r'%s, $\Delta$PMF = %1.0f $\pm$ %1.0f kJ/mol' % (name,dPMF,ddPMF))
    else:
        #plt.plot(dist,E_zero,linewidth=linewidth,color=color,label=r'%s, $\Delta$PMF = %1.0f kJ/mol' % (protein_name,dPMF))
        plt.plot(dist,E_zero,linewidth=linewidth,color=color,label=r'%s, $\Delta$PMF = %1.0f kJ/mol' % (name,dPMF))

xlim,ylim = [3.0,8.0],[-175,10]
plt.plot(xlim,[0,0],linestyle='--',linewidth=linewidth,color='grey',zorder=0)
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel('Distance [nm]',fontsize=fontsize)
plt.ylabel('PMF [kJ/mol]',fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.legend(frameon=False,loc='lower right',fontsize=fontsize)
plt.tight_layout()
plt.savefig('PMF_B.pdf',format='pdf')
plt.show()
