import numpy as np
import matplotlib.pyplot as plt

C2=1 # plot C2 and Pt dist or not?

paths = ['4/PCPSPIP2_13/CG2AT/FINAL','5/PCPSPIP2_13/CG2AT/FINAL']
colors = ['darkblue','firebrick']
names = ['PTEN','SHIP2']

if C2:
    fig,ax = plt.subplots(4,1)
else:
    fig,ax = plt.subplots(2,1)
#fig,ax = plt.subplots(2,1)

xlim = [0,600]
for path,color,name in zip(paths,colors,names):
    t,rmsd = np.genfromtxt('%s/rmsd.xvg' % path,skip_header=18+1,usecols=[0,1],unpack=True)
    t,dist = np.genfromtxt('%s/dist.xvg' % path,skip_header=17+1,usecols=[0,1],unpack=True)
    if C2:
        t,dist_C2 = np.genfromtxt('%s/dist_C2.xvg' % path,skip_header=17+1,usecols=[0,1],unpack=True)
        t,dist_Pt = np.genfromtxt('%s/dist_Pt.xvg' % path,skip_header=17+1,usecols=[0,1],unpack=True)
    idx = np.where(t>100)
    rmsd_mean = np.mean(rmsd[idx])
    dist_mean = np.mean(dist[idx])

    ax[0].plot(t,rmsd,color=color,label='%s, RMSD = %1.2f nm' % (name,rmsd_mean))

    ax[1].plot(t,dist,color=color)
    ax[1].plot(xlim,[dist[0],dist[0]],linestyle='--',color=color)
    
    if C2:
        ax[2].plot(t,dist_C2,color=color)
        ax[2].plot(xlim,[dist_C2[0],dist_C2[0]],linestyle='--',color=color)
        ax[3].plot(t,dist_Pt,color=color)
        ax[3].plot(xlim,[dist_Pt[0],dist_Pt[0]],linestyle='--',color=color)

    # print
    print('\n')
    print('protein   = %s' % name)
    print('rmsd_mean = %1.2f nm' % rmsd_mean)
    print('dist_mean = %1.2f nm' % dist_mean)

ax[0].set_ylabel('RMSD [nm]')
ax[0].set_xlim(xlim)
ax[0].set_ylim(0,1.0)
ax[0].legend(frameon=False)

if C2:
    ylim=[3.0,6.0]
else:
    ylim=[3.0,5.5]

ax[1].set_ylabel('Distance [nm]')
ax[1].set_xlim(xlim)
ax[1].set_ylim(ylim)

if C2:
    ax[2].set_ylabel('C2 Distance [nm]')
    ax[2].set_xlim(xlim)
    ax[2].set_ylim(ylim)
    ax[3].set_ylabel('Pt Distance [nm]')
    ax[3].set_xlim(xlim)
    ax[3].set_ylim(ylim)
    ax[3].set_xlabel('Time [ns]')
else:
    ax[1].set_xlabel('Time [ns]')
#ax[1].set_xlabel('Time [ns]')
#plt.tight_layout()
plt.savefig('AT')
plt.show()
