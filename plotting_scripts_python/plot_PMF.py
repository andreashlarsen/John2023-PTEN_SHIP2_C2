import numpy as np
import matplotlib.pyplot as plt

BS = 1 # has gmx wham been run with bootstrap (to get errors?)

# elements columns: protein_folder, repeat, color, time_umbrella_sampling
elements = [\
#        [1,  1, 'black',       1000],\
#        [1, 11, 'grey',        1000],\
        [2, 10, 'red',         1000, 7],\
        [2,  8, 'darkred',     1000, 7],\
        [2,  3, 'pink',        1000, 7],\
        [2, 14, 'magenta',     1000, 7]\
#        [3, 23, 'forestgreen', 1000],\
#        [3,  9, 'green',       1000],\
#        [3,  8, 'lime',        1000],\
        ]

elements_FL = [\
        [4, 13, 'darkblue',    1000, 37],\
        [5, 13, 'firebrick',   2000,  7]\
        ]

#proteins = [1,1,2,2,2,2,3,3,4,4,4,4] # 1: PTEN_Ptase, 2: SHIP2_Ptase, 3: PTEN_Ptase2, 4: PTEN_FL2, SHIP2_FL
#reps = [1,11,10,8,3,14,23,9,17,8,8,13]
#colors = ['black','grey','darkred','red','pink','magenta','forestgreen','green','cyan','blue','darkblue','deepskyblue']
#skip_first = [39,39,39,39,39,39,39,39,39,39,39]
#skip_last = [7,7,7,7,7,7,7,7,7,7,7]
#times = [1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,2000,1000]

## PROTEIN = PTEN_FL2, LIPID: PCPSPIP2
# mode 1: rep  17, frame 4505, Rzz  0.89, Dist 3.76
# mode 2: rep   8, frame 4802, Rzz -0.81, Dist 3.96
# mode 3: rep  13, frame 5336, Rzz -0.99, Dist 3.61

## PROTEIN = PTEN_Ptase, LIPID: PCPSPIP2
# mode 1: rep  1, frame 2894, Rzz -0.65, Dist 3.77
# mode 2: rep 11, frame  255, Rzz  0.99, Dist 4.10
# mode 3: rep  7, frame 3022, Rzz -0.25, Dist 3.73
# mode 4: rep  0, frame 5651, Rzz -0.51, Dist 3.66

## PROTEIN = SHIP2_Ptase, LIPID: PCPSPIP2
# mode 1: rep 20, frame 7871, Rzz -0.95, Dist 4.12
# mode 2: rep 10, frame 3781, Rzz -0.51, Dist 3.78
# mode 3: rep  8, frame 6881, Rzz  0.09, Dist 3.99
# mode 4: rep  3, frame 7346, Rzz  0.95, Dist 4.46
# mode 5: rep 14, frame 3372, Rzz -0.30, Dist 4.40

## PROTEIN = PTEN_Ptase2 (with N-terminal helix), LIPID: PCPSPIP2
# mode 1: rep 23, frame 4672, Rzz  0.77, Dist 3.68
# mode 3: rep  9, frame 3875, Rzz -0.29, Dist 3.94 some error with nan in profile.xvg - look into this
linewidth = 1.5
fontsize = 12.7
for e in elements_FL:
    prot  = e[0]
    rep   = e[1]
    color = e[2]
    time  = e[3]
    sf    = 39
    sl    = e[4]

    if prot == 1:
        protein_name = 'PTEN$_\mathrm{Pt}$'
        if rep == 1:
            mode = 1
        elif rep == 11:
            mode = 2
    elif prot == 2:
        protein_name = 'SHIP2$_\mathrm{Pt}$'
        if rep == 10:
            mode = 2
        if rep == 8:
            mode = 3
        if rep == 3:
            mode = 4
        if rep == 14:
            mode = 5
    elif prot == 3:
        protein_name = 'PTEN$_\mathrm{Pt2}$'
        if rep == 23:
            mode = 1
        if rep == 9:
            mode = 3
        if rep == 8:
            mode = 4
    elif prot == 4:
        #protein_name = 'PTEN$_\mathrm{FL2}$'
        protein_name = 'PTEN'
        if rep == 17:
            mode = 1
        if rep == 8:
            mode = 2
        if rep == 13:
            mode = 3
        k=1000
    elif prot == 5:
        #protein_name = 'SHIP2$_\mathrm{FL}$'
        protein_name = 'SHIP2'
        if rep == 13:
            mode = 1
        k=2000
    else:
        print('ERROR: protein has to be 1, 2, 3, 4 or 5')
        exit()
    
    if BS:
        filename = '%d/PCPSPIP2_%d/umbrella_t%d_l5_k2000_n1/bsResult.xvg' % (prot,rep,time)
        #dist,E,dE = np.genfromtxt(filename,skip_header=18+sf,skip_footer=sl,unpack=True)
        dist,E,dE = np.genfromtxt(filename,skip_header=18,unpack=True)
    else:
        filename = '%d/PCPSPIP2_%d/umbrella_t%d_l5_k2000_n1/profile.xvg' % (prot,rep,time)
        dist,E = np.genfromtxt(filename,skip_header=17+sf,skip_footer=sl,unpack=True)    
    zero = np.mean(E[-20:])
    E_zero = E-zero
    dPMF = np.amin(E_zero)
    #plt.plot(dist,E_zero,marker='*',color=color)
    #time_us = time/1000
    #plt.plot(dist,E_zero,color=color,label=r'%s mode %d: $\Delta$PMF = %1.0f kJ/mol,time = %d $\mu$s' % (protein_name,mode,dPMF,time_us))
    #plt.plot(dist,E_zero,linewidth=linewidth,color=color,label=r'%s mode %d: $\Delta$PMF = %1.0f kJ/mol' % (protein_name,mode,dPMF))
    
    if BS:
        plt.fill_between(dist,E_zero-dE,E_zero+dE,color=color,alpha=0.5)
        #plt.plot(dist[-20:],E_zero[-20:],marker='.',color='green')
        idx_min=np.argmin(E_zero)
        #plt.plot(dist[idx],E_zero[idx],marker='.',color='green')
        dd_zero = np.sqrt(np.sum(dE[-20:]**2))/np.sqrt(20)
        dd_min  = dE[idx_min]
        ddPMF   = np.sqrt(dd_zero**2+dd_min**2)
        plt.plot(dist,E_zero,linewidth=linewidth,color=color,label=r'%s, $\Delta$PMF = %1.0f $\pm$ %1.0f kJ/mol' % (protein_name,dPMF,ddPMF))
    else:    
        plt.plot(dist,E_zero,linewidth=linewidth,color=color,label=r'%s, $\Delta$PMF = %1.0f kJ/mol' % (protein_name,dPMF))

xlim = [3.0,8.0]
ylim = [-175,10]
plt.plot(xlim,[0,0],linestyle='--',linewidth=linewidth,color='grey',zorder=0)
plt.xlabel('Distance [nm]',fontsize=fontsize)
plt.ylabel('PMF [kJ/mol]',fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(frameon=False,fontsize=fontsize)
plt.xlim(xlim)
plt.ylim(ylim)
plt.tight_layout()
plt.savefig('PMF_A.pdf',format='pdf')
plt.show()
