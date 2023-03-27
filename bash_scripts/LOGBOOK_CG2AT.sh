#!/bin/sh
#

# --------------------------------------------------------------
# 
# This bash script runs simulations with C2 domains from different proteins and the following membranes
#
# C2 + POPC:POPS:PIP2:PIP3 (80:13:5:2) Bilayer
#
# C2 + POPC:POPS:PIP2(80:15:5) Bilayer
#
# C2 + POPC:POPS (80:20) Bilayer
#
# C2 + POPC Bilayer
# 
# CG MARTINI 2.1
# 
# Andreas Haahr Larsen
#
# November 2019
# 
# --------------------------------------------------------------

# ensure modules and other settings from batchrc file
source ~/.bashrc

###################### expect user input from here ############################################

# select frame to CG
protein_folder=3   # select protein (1=PTEN_Ptase, 2=SHIP2_Ptase, 3=PTEN_Ptase2, 4=PTEN_FL2)
lip=2              # select lipid (1=PC:PS:PIP2:PIP3, 2=PC:PS:PIP2, 3=PC:PS, 4=PC)
rep=23              # select repetition
frame=4180         # select frame (index with 1)

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

# new mode 2: rep 20, frame 7607, Rzz 0.05, dist 4.07
# new mode 3: rep 17, frame 2226, Rzz 0.99, dist 3.96

## PROTEIN = PTEN_Ptase2 (with N-terminal helix), LIPID: PCPSPIP2
# mode 1: rep 23, frame 4672, Rzz  0.77, Dist 3.68
# mode 3: rep  9, frame 3875, Rzz -0.29, Dist 3.94 
# mode 4: rep  8, frame  800, Rzz  0.99, Dist 4.40

# new mode 1: rep 7, frame 5041, Rzz -0.87, dist 3.94
# new mode 2: rep 23, frame 4180, Rzz 0.29, dist 3.67
# new mode 3: rep 21, frame 2377, Rzz 0.76, dist 3.87

## activate/deactivate (1/0) modules of script
PREPARE=0         # generate index file
BM=1              # backmapping

# Number of CPUs (for parallelization), 0 for all available
NCPUs=2

########################## to here ###################################################

# define paths to mdp files
prod=/sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/production.mdp

# set parameters depending on protein group
NEUTRAL_PROTEIN=0 # by default assume that protein has some net charge
if [ $protein_folder -eq 1 ]
then
  protein_name=PTEN_Ptase
  pdb_name=1d5r_ptase
elif [ $protein_folder -eq 2 ]
then
  protein_name=SHIP2_Ptase 
  pdb_name=5okm_B_ptase_clean
elif [ $protein_folder -eq 3 ]
then
  protein_name=PTEN_Ptase2
  pdb_name=modeller/Modeller/PTEN_Ptase2
elif [ $protein_folder -eq 4 ]
then
  protein_name=PTEN_FL2
  pdb_name=modeller/Modeller/PTEN_FL2
elif [ $protein_folder -eq 5 ]
then
  protein_name=SHIP2_FL
  pdb_name=SHIP2_FL
fi

# define path to pdb file
protein=/sansom/s157/bioc1642/Desktop/prj_PTEN_SHIP/Structures/${pdb_name}.pdb # truncated

# go to protein working directory
cd $protein_folder

# loop over lipid compositions (in this case only 1)
for j in $(seq $lip $lip)
do

  # folder name
  if [ $j -eq 1 ]
  then
    Folderprefix=PCPSPIP2PIP3
  elif [ $j -eq 2 ]
  then
    Folderprefix=PCPSPIP2
  elif [ $j -eq 3 ]
  then
    Folderprefix=PCPS
  elif [ $j -eq 4 ]
  then
    Folderprefix=PC
  fi

  # loop over replicas (in this case only 1)
  for i in $(seq $rep $rep)
  do
    
    # define folder name
    folder=${Folderprefix}_$i
    
    # go to working directory
    cd $folder
    
    ############ PREPARE BACKMAPING ###########################################################
    if [ $PREPARE -eq 1 ]
    then

      # make index file
      if [ $j -eq 1 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 28 SOL
r POPC POPS POP2 POP3
name 29 LIP
a BB
r POPS
name 31 PS
r POP2
name 32 PIP2
r POP3
name 33 PIP3
q
EOF
      elif [ $j -eq 2 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 26 SOL
r POPC POPS POP2
name 27 LIP
a BB
r POPS
name 29 PS
r POP2
name 30 PIP2
q
EOF
      elif [ $j -eq 3 ]
      then 
        cat << EOF > index.input
r W WF NA+ Ion
name 24 SOL
r POPC POPS
name 25 LIP
a BB
r POPS
name 27 PS
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 1 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 16 SOL
r POPC
name 17 LIP
a BB
q
EOF
      elif [ $j -eq 4 -a $NEUTRAL_PROTEIN -eq 0 ]
      then
        cat << EOF > index.input
r W WF NA+ Ion
name 22 SOL
r POPC
name 23 LIP
a BB
q
EOF
      fi
      gmx make_ndx -f ${protein_name}_Bilayer.gro -quiet < index.input
      rm index.input

      # see overview of index file 
      # gmx make_ndx -f ${protein_name}_Bilayer.gro -n index.ndx -quiet
  
    # end PREPARE if statement
    fi

    ############ BACKMAPPING ##################################################################
    if [ $BM -eq 1 ]
    then    

      # generate key frame index file
      cat << EOF > frame_index.ndx
[ frames ]

$frame

EOF
      # extract and center key frame from trajectory
      #echo 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -o key_frame.gro -fr frame_index.ndx -quiet
      echo 1 0 | gmx trjconv -f md.xtc -s md.tpr -n index.ndx -o key_frame.gro -fr frame_index.ndx -center -pbc mol -quiet
	  
      # generate tpr file from key frame
      gmx grompp -f $prod -c key_frame.gro -r key_frame.gro -p topol.top -o key_frame.tpr -quiet -n index.ndx

      # backmapping
      if [ $NCPUs -eq 0 ]
      then
        /sansom/s157/bioc1642/Desktop/Scripts/cg2at_owen/cg2at/cg2at -c key_frame.tpr -a $protein -ff charmm36-mar2019-updated -fg martini_2-2_charmm36 -w tip3p -loc CG2AT
      else
        /sansom/s157/bioc1642/Desktop/Scripts/cg2at_owen/cg2at/cg2at -c key_frame.tpr -a $protein -ff charmm36-mar2019-updated -fg martini_2-2_charmm36 -w tip3p -ncpus $NCPUs -loc CG2AT
      fi
    fi
    
    ############ FINISH #####################################################################

    # clean up
    rm \#*

    # navigate back to protein working directory
    cd ..

  # end loop over replicas
  done

# end loop over lipid compositions
done

# navigate back to parent directory
cd ..




