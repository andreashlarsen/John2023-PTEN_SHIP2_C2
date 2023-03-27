#!/bin/sh
#
# --------------------------------------------------------------
# 
# PROTEIN + POPC:POPS:PIP2 (80:15:5) Bilayer
#
# Umbrella sampling 
#
# CG MARTINI 2.1
# 
# Andreas Haahr Larsen
#
# February 2020
# 
# --------------------------------------------------------------

## ensure modules are loades etc from batch file
source ~/.bashrc

## set parameters

# protein folder
folder_prot=5 # (1 = PTEN_Ptase, 2 = SHIP2_Ptase, 3 = PTEN_Ptase2 (with N-term helix), 4 = PTEN_FL2 (with N-term helix), 5 = SHIP2_FL (folder 28 in prj_C2))
lip=2 # 2:PCPSPIP2, 3:PCPS, 4:PC
rep=13 
key_frame=5336

## PROTEIN = PTEN_FL2, LIPID: PCPSPIP2
# mode 1: rep  17, frame 4505, Rzz  0.89, Dist 3.76
# mode 2: rep   8, frame 4802, Rzz -0.81, Dist 3.96
# mode 3: rep  13, frame 5336, Rzz -0.99, Dist 3.61 this mode is almost identical to PTEN_FL (no helix), see slide PTEN_SHIP2_version4

## PROTEIN = PTEN_Ptase, LIPID: PCPSPIP2
# mode 1: rep  1, frame 2894, Rzz -0.65, Dist 3.77
# mode 2: rep 11, frame  255, Rzz  0.99, Dist 4.10
# mode 3: rep  7, frame 3022, Rzz -0.25, Dist 3.73
# mode 4: rep  0, frame 5651, Rzz -0.51, Dist 3.66

## PROTEIN = SHIP2_Ptase, LIPID: PCPSPIP2
# mode 1: rep  20, frame 7871, Rzz -0.95, Dist 4.12
# mode 2: rep  10, frame 3781, Rzz -0.51, Dist 3.78
# mode 3: rep   8, frame 6881, Rzz  0.09, Dist 3.99
# mode 4: rep   3, frame 7346, Rzz  0.95, Dist 4.46
# mode 5: rep  14, frame 3372, Rzz -0.30, Dist 4.40

## PROTEIN = PTEN_Ptase2 (with N-terminal helix), LIPID: PCPSPIP2
# mode 1: rep 23, frame 4672, Rzz  0.77, Dist 3.68
# mode 3: rep  9, frame 3875, Rzz -0.29, Dist 3.94
# mode 4: rep  8, frame  800, Rzz  0.99, Dist 4.40

## PROTEIN = SHIP2_FL, LIPID: PCPSPIP2
# mode 1: rep 13, frame5571, Rzz 0.99, Dist 4.26

# activate/deactivate modules of script
PREPARE=0
PULL=0
PUSH=0
GROMPP=0
MDRUN=0
AXON=0
DOWNLOAD=0
ANALYSIS=1

# define paths to mdp files
umb=/sansom/s157/bioc1642/Desktop/Scripts/umbrella/umbrella.mdp
pull=/sansom/s157/bioc1642/Desktop/Scripts/umbrella/pull.mdp

# define paths to scripts 
py3=/sansom/s157/bioc1642/anaconda3/bin/python3.7
extract_frames=/sansom/s157/bioc1642/Desktop/Scripts/umbrella/extract_frames.py

# parallelisation with mdrun
pinoffset=18
noCPUs=6
noGPUs=1

# cpu for cpu only, auto for using gpu
nb=auto

# parameters for push and pull and umbrella
if [ $folder_prot -eq 5 ] 
then
time=2000 # simulation time in umbrella sampling
elif [ $folder_prot -eq 4 ]
then
time=1000 # simulation time in umbrella sampling
else
time=2000
echo ----------------------------
echo time set to default: 2000 ns
echo ----------------------------
fi

step_size=5
# 50 = 0.5 nm
# 20 = 0.2 nm
# 10 = 0.1 nm
#  5 = 0.05 nm
#  2 = 0.02 nm
#  1 = 0.01 nm
umb_rep=1
spring_const=2000
dt=35 # timestep, in fs, in umbrella.mdp
nsteps=$((time*1000*1000/dt))
echo nsteps in umbrella is $nsteps

# set max dist in extract_frames (for umbrella sampling)
if [ $folder_prot -eq 4 ] || [ $folder_prot -eq 5 ]
then
  max_dist=8.5
else
  max_dist=7.0
fi
echo max_dist in extract_frames.py is $max_dist

# folder name
if [ $lip -eq 1 ]
then
  Folderprefix=PCPSPIP2PIP3
elif [ $lip -eq 2 ]
then
  Folderprefix=PCPSPIP2
elif [ $lip -eq 3 ]
then
  Folderprefix=PCPS
elif [ $lip -eq 4 ]
then
  Folderprefix=PC
fi
echo Folderprefix = lipid composition is $Folderprefix
    
############ GENERATE FOLDERS ETC #######################################################
    
## define folder name
folder=${folder_prot}/${Folderprefix}_$rep
cd $folder

## go to folder (directory)
dir=umbrella_t${time}_l${step_size}_k${spring_const}_n${umb_rep}
mkdir -p $dir
cd $dir
    
############ PREPARE TO RUN SIMULATION ###################################################

if [  $PREPARE -eq 1 ]
then

echo --------------------------------------------
echo PREPARE to run sim
echo --------------------------------------------
	
## extract key frame from trajectory
cat << EOF > frame_index.ndx
[ frames ]

$key_frame

EOF
echo 0 | gmx trjconv -f ../md.xtc -s ../md.tpr  -o key_frame.gro -fr frame_index.ndx -quiet
       
## make index file
if [ $lip -eq 2 ]
then
  cat << EOF > index.input
r W WF ION
name 19 SOL
r POPC POPS POP2
name 20 LIP
q
EOF
elif [ $lip -eq 3 ]
then
  cat << EOF > index.input
r W WF ION
name 18 SOL
r POPC POPS
name 19 LIP
q
EOF
elif [ $lip -eq 4 ]
then
  cat << EOF > index.input
r W WF ION
name 17 SOL
r POPC
name 18 LIP
q
EOF
fi
gmx make_ndx -f ../md.gro -quiet < index.input
rm index.input

fi 

############ PULL SIMULATION ###################################################

if [ $PULL -eq 1 ]
then

  echo --------------------------------------------
  echo make PULL sim
  echo --------------------------------------------

  ## make pull simulation
  cp $pull pull_slow.mdp
  sed -i -e "s/pull_coord1_rate        = 0.0001/pull_coord1_rate        = 0.00002/g" pull_slow.mdp
  ## restrain POPC and POPS movement (still restrain POP2)
  sed -i -e "s/define                   = -DPOSRES_POP2/define                   = -DPOSRES_POP2 -DPOSRES_POPC -DPOSRES_POPS/g" pull_slow.mdp
  gmx grompp -f pull_slow.mdp -c key_frame.gro -r key_frame.gro -p ../topol.top -n index.ndx -o pull.tpr -quiet -maxwarn 1
  gmx mdrun -deffnm pull -v -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet -nsteps 1000 # very short run, to get pull.gro
  gmx mdrun -deffnm pull -v -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet # it will run until it crashes (after roughly 20 min)
    
  ## center and fix PBC for vizualization of trajectory
  echo 1 0 | gmx trjconv -f pull.xtc -s pull.tpr -dt 1000 -pbc mol -center -o pull_dt1000.xtc -quiet
  cat << EOF > pull.pml
load pull.gro
load_traj pull_dt1000.xtc
remove resname W
remove resname WF
remove resname ION
extract LIP, resname POPC resname POPS resname POP2
show spheres
spectrum count, rainbow, pull
show cell
set opaque_background, off
set orthoscopic, on
EOF
  #pymol pull.pml
  
  ## check position vs time and force vs time
  #xmgrace pull_pullx.xvg
  #xmgrace pull_pullf.xvg  
   
  ## extract frames from steered MD trajectory 
  cp $extract_frames .
  sed -i -e "s/step_size = 5  /step_size = ${step_size}  /g" extract_frames.py # change step size
  echo changing max dist in extract_frames.py
  sed -i -e "s/max_dist = 7.0/max_dist = ${max_dist}/g" extract_frames.py # change max dist
  $py3 extract_frames.py # generates extract_frames.ndx
  echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -n  -fr extract_frames.ndx -o extract_frames.xtc -quiet
  ## check number of extracted frames
  #gmx check -f extract_frames.xtc -quiet

fi

############ PUSH SIMULATION ###################################################

if [ $PUSH -eq 1 ]
then

  echo --------------------------------------------
  echo make PUSH sim
  echo --------------------------------------------      

  ## make push simulation
  cp $pull push_slow.mdp
  sed -i -e "s/pull_coord1_rate        = 0.0001/pull_coord1_rate        = -0.00002/g" push_slow.mdp
  # no restrain on lipids
  sed -i -e "s/define                   = -DPOSRES_POP2/;define                   = -DPOSRES_POP2 -DPOSRES_POPC -DPOSRES_POPS/g" push_slow.mdp
  gmx grompp -f push_slow.mdp -c key_frame.gro -r key_frame.gro -p ../topol.top -n index.ndx -o push.tpr -quiet -maxwarn 1
  
  if [ $folder_prot -eq 4 ]
  then
    gmx mdrun -deffnm push -v -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet -nsteps 4000000
    # nsteps of 4000000 corresponds to 150  ns, and 3.0 nm with rate of 0.00002 nm/ps and dt of 30 fs.
    echo push protein com 2.4 into membrane
  else
    gmx mdrun -deffnm push -v -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet -nsteps 2500000
    # nsteps of 2500000 corresponds to 75  ns, and 1.5 nm with rate of 0.00002 nm/ps and dt of 30 fs.
    echo push protein com 1.5 into membrane
  fi
  
  ## center and fix PBC for vizualization of trajectory
  echo 1 0 | gmx trjconv -f push.xtc -s push.tpr -dt 1000 -pbc mol -center -o push_dt1000.xtc -quiet
cat << EOF > push.pml
load push.gro
load_traj push_dt1000.xtc
remove resname W
remove resname WF
remove resname ION
extract LIP, resname POPC resname POPS resname POP2
show spheres
spectrum count, rainbow, push
show cell
set opaque_background, off
set orthoscopic, on
EOF
  #pymol push.pml

  ## check position vs time and force vs time
  #xmgrace push_pullx.xvg
  #xmgrace push_pullf.xvg

  ## extract frames from steered MD trajectory
  cp $extract_frames extract_frames_push.py
  sed -i -e "s/step_size = 5  /step_size = ${step_size}  /g" extract_frames_push.py # change step size
  sed -i -e "s/file_in = 'pull_pullx.xvg'/file_in = 'push_pullx.xvg'/g" extract_frames_push.py
  sed -i -e "s/file_out = 'extract_frames.ndx'/file_out = 'extract_frames_push.ndx'/g" extract_frames_push.py  
  sed -i -e "s/min_dist = dist\[0\]/min_dist = dist[-1]/g" extract_frames_push.py
  sed -i -e "s/max_dist = 7.0/max_dist = dist[0]/g" extract_frames_push.py
  sed -i -e "s/np.where(dist>d)/np.where(dist<d)/g" extract_frames_push.py

  $py3 extract_frames_push.py # generates extract_frames.ndx
  echo 0 | gmx trjconv -f push.xtc -s push.tpr -n  -fr extract_frames_push.ndx -o extract_frames_push.xtc -quiet

  ## check number of extracted frames
  #gmx check -f extract_frames.xtc -quiet
fi


############ COUNT NUMBER OF FRAMES ###################################################

if [ $GROMPP -eq 1 ] || [ $MDRUN -eq 1 ] || [ $AXON -eq 1 ]
then

  ## get number of frames (3 lines in file are not frames numbers)
  Nframes_pull=-3; while read -r LINE; do (( Nframes_pull++ )); done < extract_frames.ndx
  echo $Nframes_pull
  Nframes_push=-3; while read -r LINE; do (( Nframes_push++ )); done < extract_frames_push.ndx
  echo $Nframes_push    
  Nframes=$((Nframes_push+Nframes_pull))
  echo $Nframes
fi

############ GROMPP (FOR UMBRELLA WINDOW SIMULATIONS) ###################################################

if [ $GROMPP -eq 1 ]
then

  echo --------------------------------------------
  echo GROMPP
  echo --------------------------------------------

  ## change production time and spring constant of umbrella run
  cp $umb umbrella.mdp
  sed -i -e "s/pull_coord1_k           = 1000/pull_coord1_k           = ${spring_const}/g" umbrella.mdp
  sed -i -e "s/nsteps               = 30000000/nsteps               = ${nsteps}/g" umbrella.mdp    
  sed -i -e "s/1050 ns/2000 ns/g" umbrella.mdp    
  ## remove restraine on lipids in production run
  sed -i -e "s/define                   = -DPOSRES_POP2/;define                   = -DPOSRES_POP2/g" umbrella.mdp

  #rm tpr_files_umbrella.dat
  #rm pullf_files_umbrella.dat
  #for frame in $(seq 89 $Nframes)
  for frame in $(seq 1 $Nframes)
  do
    echo --------------------------------
    echo GROMPPing frame $frame of $Nframes
    echo --------------------------------

    if [ $frame -le $Nframes_push ]
    then 
	cat << EOF > frames_step.ndx
 [ frames ]

$frame

EOF
	echo 0 | gmx trjconv -f extract_frames_push.xtc -s push.tpr -fr frames_step.ndx -o step_${frame}.gro -quiet
    else
	cat << EOF > frames_step.ndx
 [ frames ]

$((frame-Nframes_push))

EOF
	echo 0 | gmx trjconv -f extract_frames.xtc -s pull.tpr -fr frames_step.ndx -o step_${frame}.gro -quiet
    fi

    ## grompp umbrella production runs
    gmx grompp -f umbrella.mdp -c step_${frame}.gro -r step_${frame}.gro -n index.ndx -p ../topol.top -o umbrella_step${frame}.tpr -maxwarn 1 -quiet
    	
    ## add lines to summary files for WHAM
    echo "umbrella_step${frame}.tpr" >> tpr_files_umbrella.dat
    echo "umbrella_step${frame}_pullf.xvg" >> pullf_files_umbrella.dat
    
  done
  
  # cp summary files - ready for truncation
  cp tpr_files_umbrella.dat tpr_files_umbrella_trunc.dat
  cp pullf_files_umbrella.dat pullf_files_umbrella_trunc.dat
fi	
    
############ RUN UMBRELLA WINDOWS ##############################################################

if [ $MDRUN -eq 1 ]
then

  echo --------------------------------------------
  echo make UMBRELLA sims
  echo --------------------------------------------

  ## mdrun
  # performance, 4 cpu, 1 gpu: 9 us/day, so 1 run should not take more than about 3h
  for frame in $(seq 1 $Nframes)
  #for frame in $(seq 89 103)
  #for frame in $(seq 104 $Nframes)
  do
    echo --------------------------------------------
    echo UMBRELLA frame $frame out of $Nframes
    echo --------------------------------------------
    # -maxh added to continue to next step after 5h if simulation crashes, which happens sometimes 
    gmx mdrun -deffnm umbrella_step${frame} -pin on -ntomp $noCPUs -ntmpi $noGPUs -pinoffset $pinoffset -quiet -nb $nb -maxh 5.0 
  done
fi

if [ $AXON -eq 1 ]
then
  ## scp tpr files to Axon
  axon_dir=axon_umbrella_prot${folder_prot}_lip${Folderprefix}_rep${rep}_t${time}_l${step_size}_k${spring_const}_n${umb_rep}
  mkdir -p $axon_dir
  cp umbrella_step*.tpr $axon_dir

  ## make sbatch file and copy to Axon
  cat << EOF > Flow_axon.sh
#!/bin/bash 

## Define if you are in testing mode (on short queue) or running

## Set the job name.
#SBATCH --job-name=pro${folder_prot}_lip${lip}

## Set the number of nodes and cores. This shouldn't need changing.
## The maximum number of nodes which can be used on axon is 1.
#SBATCH --nodes=1

## Set the number of tasks on each node, this is usually the number of program executions done.
## For example, if you are running an mpi run, this would reflect the number passed to the -np flag.
#SBATCH --ntasks-per-node=1

## Set the number of cores per tasks. This usually reflects the number of threads (often openmp) that
## are being assigned per task. Benchmarking has shown that setting tasks to 1 and cpus-per-task to
## 16 for atomistic simulations and 8 for CG is a good starting point for getting maximum efficiency.
#SBATCH --cpus-per-task=16

## Set the number of GPUs to be used. In most cases this will be set to 1.
#SBATCH --gres=gpu:1

## IMPORTANT: set GPU binding, otherwise your jobs will clash with other users'
#SBATCH --gres-flags=enforce-binding

## Select the queues you will be running on (sansom: gpu-sansom,gpu-sansom2 biggin: gpu-biggin,gpu-biggin2) 
##SBATCH -p gpu-sm-short
#SBATCH -p gpu-sansom,gpu-sansom2

## Select the max amount of time this job will run (48h for gpu-sansom, 3h for shor)
##SBATCH --time=3:00:00
#SBATCH --time=48:00:00

## Set an array of jobs you want to run
#SBATCH --array=1-$Nframes%12

source /etc/profile.d/modules.sh

module purge
module load apps/gromacs/2018.6-plumed_2.4.4-GPU-KEPLER

## Note if running an energy minimisation, CG or using energy groups, you need to add the -ntmpi 1 for gromacs 2019 and above
gmx mdrun -deffnm umbrella_step\${SLURM_ARRAY_TASK_ID} -v -ntomp \${SLURM_CPUS_PER_TASK} -ntmpi 1 -quiet
EOF
  cp Flow_axon.sh $axon_dir

  scp -r $axon_dir axon:/home/bioc1642/
  # remove folder after scp
  rm -r $axon_dir
fi
   
if [ $DOWNLOAD -eq 1 ]
then
  axon_dir=axon_umbrella_prot${folder_prot}_lip${Folderprefix}_rep${rep}_t${time}_l${step_size}_k${spring_const}_n${umb_rep}
  #scp -r axon:/home/bioc1642/$axon_dir .
  #cp $axon_dir/umbrella_step*_pullf.xvg .
  scp -r axon:/home/bioc1642/$axon_dir/umbrella_step*_pullf.xvg .

  ## vizualization
  name=umbrella_step100
  # download 
  scp -r axon:/home/bioc1642/$axon_dir/$name.gro .
  scp -r axon:/home/bioc1642/$axon_dir/$name.xtc .
  # center and fix PBC
  echo 1 0 | gmx trjconv -f $name.xtc -s $name.tpr -dt 1000 -pbc mol -center -o ${name}_dt1000.xtc -quiet
  rm umbrella_step20.xtc
  # pymol script 
  cat << EOF > umbrella.pml
load $name.gro
load_traj ${name}_dt1000.xtc
remove resname W
remove resname WF
remove resname ION
extract LIP, resname POPC resname POPS resname POP2
show spheres
spectrum count, rainbow, $name
show cell
set opaque_background, off
set orthoscopic, on
EOF
  #pymol umbrella.pml
fi

############ ANALYSIS ###################################################################
if [ $ANALYSIS -eq 1 ]
  then

  ## run WHAM (weighted histogram analysis method)
  #skip=$(( time*1000/5*0 )) # skip first fifth of frames (equilibration) for WHAM
  #gmx wham -tol $tol -it tpr_files_umbrella_trunc.dat -if pullf_files_umbrella_trunc.dat -hist -o -bsres -temp 323 -quiet -b $skip -e 80000e3 #-nBootstrap 200 # outcomment "-nBootstrap 200" to get uncertainties 

  if [ $folder_prot -eq 5 ]
  then
    taskset -c 0-1 gmx wham -it tpr_files_umbrella_trunc.dat -if pullf_files_umbrella_trunc.dat -hist -temp 323 -quiet -b 100 -e 10050000 -nBootstrap 200 -min 3.3 -max 8 # only up to 1050 ns - to make it the same as the partial_PMFs
  elif [ $folder_prot -eq 4 ]
  then
    taskset -c 0-1 gmx wham -it tpr_files_umbrella_trunc.dat -if pullf_files_umbrella_trunc.dat -hist -temp 323 -quiet -b 100 -e 10050000 -nBootstrap 200 -min 3.3 -max 8
  else
    echo ----------------------------
    echo not running WHAM
    echo ----------------------------
  fi

  #taskset -c 0-1 gmx wham -it tpr_files_umbrella_trunc.dat -if pullf_files_umbrella_trunc.dat -hist -temp 323 -quiet -b 100000 -e 10050000 -nBootstrap 200 -min 3 -max 7.5 	# only up to 1050 ns - to make it the same as the partial_PMFs 
  #taskset -c 0-1 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 323 -quiet -b 100000 -e 10050000 -nBootstrap 200 -min 3 -max 7.5 	# only up to 1050 ns - to make it the same as the partial_PMFs 
  #gmx wham -tol $tol -it tpr_files_umbrella_trunc.dat -if pullf_files_umbrella_trunc.dat -hist -o -bsres -temp 323 -quiet -b $skip #-nBootstrap 200 # outcomment "-nBootstrap 200" to get uncertainties  
      
  # check PMF curve
  xmgrace -nxy histo.xvg
  #xmgrace profile.xvg
  xmgrace bsProfs.xvg
  xmgrace bsResult.xvg
fi

############ FINISH #####################################################################

## clean up
rm \#*

## navigate back to parent directory
cd ../../..

