#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Rzz against distance
Andreas Haahr Larsen
"""

# import libraries
import numpy as np
import pylab
import matplotlib as mpl
from matplotlib import pyplot as plt

#fontsize
plt.rcParams.update({'font.size': 12  })
    
# protein folder (1=PTEN_Ptase, 2=SHIP2_Ptase, 3=PTEN_Ptase2, 4=PTEN_FL2)
Folder = 3

# plot only selected trajectories (or a single trajectory)
SELECTED_TRAJS = 1
trajs = [3,21]

## PROTEIN = SHIP2_Ptase, LIPID: PCPSPIP2
# trajs = [0] # mode 1
# trajs = [3,5,19,20] # mode 2
# trajs = [2,4,6,7,8,9,12,13,14,15,16,17,18,22,24] # mode 3 (9 and 17 has Rzz very close to 1)

# trajs = [10,23] # mode 4, (between mode 2 and 3
# trajs = [11] # mode 2 and 4
# trajs = [1] # mode 1->2
# trajs = [21] # mode 3 and 4

# selected frames:
# mode 2. rep 20, frame 7607
# mode 3: rep 17, frame 2226

## PROTEIN = PTEN_Ptase2 (with helix), LIPID: PCPSPIP2
# trajs = [7,9] # mode 1 
# trajs = [0,2,8,10,11,12,14,15,16,18,20,22,23,24] # mode 2 
# trajs = [3,21] # mode 3

# trajs = [1,6] # mode 1 and 3
# trajs = [5,13,17] # mode 2 and 3
# trajs = [4,19] # between mode 2 and 3

# selected frames:
# mode 1: rep 7, frame 5041, Rzz -0.87, dist 3.94 
# mode 2: rep 23, frame 4180, Rzz 0.29, dist 3.67 
# mode 3: rep 21, frame 2377, Rzz 0.76, dist 3.87

# bins
bins = 100

# find frame corresponding to point on Dist vs. Rzz plot
aTol = 2 # increase tolerance area by a factor 
RzzInt = 0.76
DistInt = 3.87
RzzTol = 0.002*aTol
DistTol = 0.002*aTol

# include only last frames 286 frames (out of 5734), i.e. 100 ns out of 2 us
LAST = 0
last_frames = 500

# choose to save figures or not
SAVE = 1
# save high resolution (and without title and grid)
HIGH_RES = 0
# save high res with color bar as Rzz_color_bar.eps
HIGH_RES_COLOR_BAR = 0 

# plot data from all repeats
PLOT_ALL = 0

# lipid compositon, color, number of repetitions
Reps = [24,0,0] # PCPSPIP2, PCPS, PC

if Folder == 1:
    Protein = "PTEN_Ptase"
elif Folder == 2:
    Protein = "SHIP2_Ptase"
elif Folder == 3:
    Protein = "PTEN_Ptase2"
elif Folder == 4:
    Protein = "PTEN_FL2"

print('\n##################################################\nProtein: %s' % Protein)

Elements = [\
            ["PCPSPIP2","red",Reps[0],"PC:PS:PIP$_2$"],\
            ["PCPS","blue",Reps[1],"PC:PS"],\
            ["PC","green",Reps[2],"PC"],\
            ] 
for Element in Elements:
    Lip = Element[0]
    Color = Element[1]
    Repeats = Element[2]
    Lip_name = Element[3]
    if Repeats > 0:
        RzzAll = []
        DistAll = []
        print('\n##################################################\nLipid: %s' % Lip)
        if SELECTED_TRAJS == 0:
            trajs = np.arange(Repeats)
        for i in trajs:
            # print('i = %d' % i)
            # Rzz
            File = str(Folder) + "/" + Lip + "_" + str(i) + "/Rzz.xvg"
            Time,Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz = np.genfromtxt(File,skip_header=32,usecols=[0,1,2,3,4,5,6,7,8,9],unpack=True)
            Time /= 1000 # convert from ps to ns
            
            if LAST:
                Time = Time[-last_frames:]
                Rzz = Rzz[-last_frames:]
                
            if PLOT_ALL:
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Time,Rzz,label="Rzz",color="r")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Time [ns]")
                plot1.set_ylabel("Rzz [nm]")
                pylab.suptitle("Rzz")
                pylab.show()
            
            RzzAll = np.append(RzzAll,Rzz)
        
            # Dist
            File = str(Folder) + "/" + Lip + "_" + str(i) + "/dist_com_cent.xvg"
            Time,Dist = np.genfromtxt(File,skip_header=17,usecols=[0,1],unpack=True)
            
            if LAST:
                Time = Time[-last_frames:]
                Dist = Dist[-last_frames:]
            
            # find structure of interest
            Rzz_filter_index = np.where((Rzz > RzzInt-RzzTol) & (Rzz < RzzInt + RzzTol))
            Dist_filter_index = np.where((Dist > DistInt-DistTol) & (Dist < DistInt + DistTol))
            Common_index = np.intersect1d(Rzz_filter_index,Dist_filter_index)
            l_Common = len(Common_index)
            
            if l_Common:
                for j in range(l_Common):
                    print('Structure of interest: lip %s, repeat %d, frame %d' % (Lip,i,Common_index[j]))
            # for debubbing:
            #else:    
                #len_Rzz = len(Rzz_filter_index[0])
                #len_Dist = len(Dist_filter_index[0])
                #print('left Rzz values = %d, left Dist values = %d, left common values = %d' % (len_Rzz,len_Dist,l_Common) )
            
            if PLOT_ALL:  
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Time,Dist,label="Dist",color="g")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Time [ns]")
                plot1.set_ylabel("Dist [nm]")
                pylab.suptitle("Dist")
                pylab.show()
                
            DistAll = np.append(DistAll,Dist)
            
            # Together
            if PLOT_ALL:
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Rzz,Dist,label="Dist vs Rzz",color="k",marker=".",linestyle="none")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Rzz")
                plot1.set_ylabel("Dist")
                pylab.suptitle("Dist vs Rzz")
                pylab.show()
            
        # Histogram Rzz vs Dist
        hist=plt.hist2d(RzzAll,DistAll,bins=bins,norm=mpl.colors.LogNorm())
        
        
        plt.clim(1,1000)
        
            
        plt.ylim((2.5,7.5))
        plt.xlim((-1.0,1.0))
        plt.tight_layout()
        if HIGH_RES == 0:
            #plt.title("C2 from " + Protein + " and " + Lip + " membrane")
            plt.minorticks_on()
            plt.grid()
            plt.grid(which='minor')
            plt.xlabel('$R_{zz}$')
            plt.ylabel('Distance [nm]')
            plt.colorbar()
        else:
            if HIGH_RES_COLOR_BAR:
                plt.colorbar()
            else:
                plt.title(r'%s' % Lip_name)
                if Folder == 15:
                    plt.xlabel('$R_{zz}$')    
                if Lip == 'PC':
                    plt.ylabel('Distance [nm]')
#            if Lip == 'PCPSPIP2':
#                plt.colorbar()    
        plt.tight_layout()     
#        if SAVE and SELECTED_TRAJS == 0:
#            if LAST:
#                plt.savefig('../../Seafile/C2_Manus/Rzz_%s_%s_last.png' % (Protein,Lip))
#            else:
#                if HIGH_RES:
#                    if HIGH_RES_COLOR_BAR:
#                        plt.savefig('Rzz_color_bar.eps', format='eps')
#                    else:
#                        plt.savefig('Rzz_%s_%s.eps' % (Protein,Lip), format='eps')
#                    
#                else:
#                    plt.savefig('Rzz_%s_%s.png' % (Protein,Lip))

        # find bin with highest occupancy
        maxbin=np.amax(hist[0])
        indc = np.where(maxbin==hist[0])
        if len(indc) > 1:
            indxRzz = indc[0][0]
            indxDist = indc[1][0]
        else:
            indxRzz = indc[0]
            indxDist = indc[1]
        Rzz_indx_min = hist[1][indxRzz]
        Rzz_indx_max = hist[1][indxRzz+1]
        Dist_indx_min = hist[2][indxDist]
        Dist_indx_max = hist[2][indxDist+1]
        Rzz_indx_mean = (Rzz_indx_min + Rzz_indx_max)/2
        Dist_indx_mean = (Dist_indx_min + Dist_indx_max)/2
        
        print('Bin with highest density: Rzz = %1.2f, Dist = %1.2f' % (Rzz_indx_mean,Dist_indx_mean))
        
        plt.show()
        if SAVE:
            plt.savefig('Rzz_%s_%s.png' % (Protein,Lip))

