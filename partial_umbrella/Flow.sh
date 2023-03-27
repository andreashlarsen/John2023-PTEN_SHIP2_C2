
## select
ITP=0
SEPARATE=0
SOLVATE=0
PULL=0
GROMPP=0
RUN_LOCAL=0
DOWNLOAD=0
ANALYZE=1

## performance
pinoffset=0
noCPUs=4

## protein
#protein=PTEN_Pt
#protein=PTEN_C2 
#protein=SHIP2_Pt 
protein=SHIP2_C2
mkdir -p ${protein}/umbrella

## axon dir
axon_dir=axon_umb_$protein

######################################################################################################
### NEW Protein.itp FILE FOR PTEN_C2 #################################################################
######################################################################################################

if [ $ITP -eq 1 ]
then
  path=/sansom/s157/bioc1642/Desktop/prj_PTEN_SHIP/Structures/modeller/Modeller/PTEN_C2/PTEN_C2.pdb
  martinize=/sansom/s157/bioc1642/Desktop/Scripts/martinize_GROMACS_2018_plumed.py # edited line 1851 to have "/gromacs/top/" instead of "/top/"
  python $martinize -f $path -x PTEN_C2_CG.pdb -o topol_PTEN_C2.top -v -ff martini22 -elastic -dssp dssp
  rm PTEN_C2_CG.pdb chain_.ssd topol_PTEN_C2.top
  mv Protein.itp Protein_PTEN_C2.itp
fi

######################################################################################################
### SEPARATE #########################################################################################
######################################################################################################

if [ $SEPARATE -eq 1 ]
then
## import FL CG structures
cp ../4/PCPSPIP2_13/key_frame.gro PTEN_FL_mode1.gro
cp ../4/PCPSPIP2_13/key_frame.tpr PTEN_FL_mode1.tpr
cp ../4/PCPSPIP2_13/topol.top topol_PTEN.top

cp ../../prj_C2/28/PCPSPIP2_13/key_frame.gro SHIP2_FL_mode1.gro
cp ../../prj_C2/28/PCPSPIP2_13/key_frame.tpr SHIP2_FL_mode1.tpr
cp ../../prj_C2/28/PCPSPIP2_13/topol.top topol_SHIP2.top

## import Ptase and C2 topology files
cp ../3/PCPSPIP2_0/Protein.itp Protein_PTEN_Pt.itp
i#cp ../../prj_C2/23/PCPSPIP2_0/Protein.itp Protein_PTEN_C2.itp
cp ../2/PCPSPIP2_0/Protein_B.itp Protein_SHIP2_Pt.itp
cp ../../prj_C2/25/PCPSPIP2_0/Protein.itp Protein_SHIP2_C2.itp

## Indices of C2 and Ptase domains

# PTEN_Ptase2: 1-186
# PTEN_linker: 187-190
# PTEN_C2: 191-352 

# SHIP2_Ptase: 1-314
# SHIP2_linker: 315-330
# SHIP2_C2: 332-457
# SHIP2_Cterm: 458-461

## extract lip+PTEN_Pt or lip+PTEN_C2
cat << EOF > index_input
ri 1-186
name 19 Ptase
ri 191-352
name 20 C2
13 | 14 | 15
name 21 LIP
21 | 19
name 22 LIP_Ptase
21 | 20
name 23 LIP_C2
q
EOF
gmx make_ndx -f PTEN_FL_mode1.gro -quiet < index_input
echo 22 | gmx trjconv -f PTEN_FL_mode1.gro -s PTEN_FL_mode1.tpr -n index -o PTEN_Pt_lip.gro -quiet
echo 23 | gmx trjconv -f PTEN_FL_mode1.gro -s PTEN_FL_mode1.tpr -n index -o PTEN_C2_lip.gro -quiet

## extract lip+SHIP2_Pt or lip+SHIP2_C2
cat << EOF > index_input
ri 3-314
name 19 Ptase
ri 332-457
name 20 C2
13 | 14 | 15
name 21 LIP
21 | 19
name 22 LIP_Ptase
21 | 20
name 23 LIP_C2
q
EOF
gmx make_ndx -f SHIP2_FL_mode1.gro -quiet < index_input
echo 22 | gmx trjconv -f SHIP2_FL_mode1.gro -s SHIP2_FL_mode1.tpr -n index -o SHIP2_Pt_lip.gro -quiet
echo 23 | gmx trjconv -f SHIP2_FL_mode1.gro -s SHIP2_FL_mode1.tpr -n index -o SHIP2_C2_lip.gro -quiet

## check index files
#gmx make_ndx -f PTEN_FL_mode1.gro -n index.ndx -quiet
#gmx make_ndx -f SHIP2_FL_mode1.gro -n index.ndx -quiet

fi

######################################################################################################
### SOLVATE AND EQUILIBRATE ##########################################################################
######################################################################################################

if [ $SOLVATE -eq 1 ] 
then

## import scripts
py2=/sansom/s157/bioc1642/anaconda2/bin/python2.7
insane=/sansom/s157/bioc1642/Desktop/Scripts/insane.py

## change dir
cd $protein
cp ../${protein}_lip.gro .
cp ../Protein_${protein}.itp .

## solvate 
if [ $protein == PTEN_Pt ] || [ $protein == PTEN_C2 ]
then
  echo Protein domain from PTEN
  $py2 $insane -f ${protein}_lip.gro -o ${protein}_solv.gro -p topol.top -x 12.22851 -y 12.22851 -z 25.08456 -sol W:9 -sol WF:1 -salt 0
  cp ../topol_PTEN.top topol_FL.top
else
  echo Protein domain from SHIP2
  $py2 $insane -f ${protein}_lip.gro -o ${protein}_solv.gro -p topol.top -x 12.09219 -y 12.09219 -z 25.68106 -sol W:9 -sol WF:1 -salt 0
  cp ../topol_SHIP2.top topol_FL.top
fi

## change topology file
head -n -3 topol_FL.top > tmp1
tail -3 topol.top > tmp2
cat tmp1 tmp2 > topol_new.top
rm tmp1 tmp2 topol.top topol_FL.top
sed -i -e "s/Protein.itp/Protein_$protein.itp/g" topol_new.top

if [ $protein == SHIP2_Pt ]
then
    sed -i -e "s/Protein        1/Protein_B        1/g" topol_new.top
fi

## minimization
cat << EOF > min.mdp
; change some constraints to bonds - necessary for PIPs
define                   = -DFLEXIBLE
; Run control
integrator               = steep
nsteps                   = 20000
; Output control
nstxout                  = 0
nstfout                  = 0
nstlog                   = 100
; Neighbour searching
cutoff-scheme            = Verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
; Electrostatics
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15    ; 2.5 (with polarizable water)
epsilon_rf               = 0
; Van der Waals
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
EOF
gmx grompp -f min.mdp -c ${protein}_solv.gro -p topol_new.top -o min.tpr -quiet
gmx mdrun -deffnm min -quiet -pin on -ntomp $noCPUs -ntmpi 1 -pinoffset $pinoffset

## make index file (with solvent)
cat << EOF > index_input
13 | 14 | 15
name 26 LIP
16 | 17 | 18
name 27 SOL
q
EOF
gmx make_ndx -f ${protein}_solv.gro -quiet < index_input

## position restraint on protein
echo 1 | gmx genrestr -f min.gro -n index.ndx -o posre.itp -quiet

## short equilibration
cat <<EOF > eq.mdp
; position restraint
define                   = -DPOSRES  
refcoord_scaling         = all 
; Run control
integrator               = md
dt                       = 0.03
nsteps                   = 333333 ; 333333*30 fs = 10 ns
; Output control
nstlog                   = 10000 ; Output frequency for energies to log file
nstenergy                = 100   ; Output frequency for energies to energy file
nstxout-compressed       = 1000  ; Output frequency for .xtc file
; Neighbour searching
cutoff-scheme            = Verlet
nstlist                  = 20     ; good performance for GPUs, default is 10
ns_type                  = grid
pbc                      = xyz
; Electrostatics
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15                ; 2.5 (with polarizable water)
; Van der Waals
vdw_type                 = cutoff
rvdw                     = 1.1
; Temperature coupling
tcoupl                   = v-rescale
tc-grps                  = Protein LIP SOL
tau_t                    = 1.0  1.0  1.0
ref_t                    = 323 323 323
Pcoupl                   = berendsen
Pcoupltype               = semiisotropic
tau_p                    = 6.0  ;       parrinello-rahman is more stable with larger tau-p, DdJ, 20130422
compressibility          = 3e-4 3e-4
ref_p                    = 1.0 1.0
; Velocity generation
gen_vel                  = yes
gen_temp                 = 323
gen_seed                 = 473529
; Bonds and constraints
constraints              = none
constraint_algorithm     = Lincs
continuation             = no
EOF
gmx grompp -f eq.mdp -c min.gro -r min.gro -p topol_new.top -o eq.tpr -n index.ndx -quiet
gmx mdrun -v -deffnm eq -quiet -pin on -ntomp $noCPUs -ntmpi 1 -pinoffset $pinoffset

cd ..

fi 

######################################################################################################
### PULL and PUSH ####################################################################################
######################################################################################################

if [ $PULL -eq 1 ]
then

## change dir
cd ${protein}/umbrella

## pull
cat <<EOF > pull.mdp
; position restraint
define                   = -DPOSRES_POP2   ; Restraint on POP2 (maybe need to restrain other lips: define                   = -DPOSRES_POP2 -DPOSRES_POPC -DPOSRES_POPS
refcoord_scaling         = all 
; Run control
integrator               = md
dt                       = 0.03
nsteps                   = 33333333 ; 33333333*30 fs = 1000 ns
; Output control
nstlog                   = 10000 ; Output frequency for energies to log file
nstenergy                = 300   ; Output frequency for energies to energy file
nstxout-compressed       = 300  ; Output frequency for .xtc file
; Neighbour searching
cutoff-scheme            = Verlet
nstlist                  = 20     ; good performance for GPUs, default is 10
ns_type                  = grid
pbc                      = xyz
; Electrostatics
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15
; Van der Waals
vdw_type                 = cutoff
rvdw                     = 1.1
; Temperature coupling
tcoupl                   = v-rescale
tc-grps                  = Protein LIP SOL
tau_t                    = 1.0  1.0  1.0
ref_t                    = 323 323 323
Pcoupl                   = berendsen
Pcoupltype               = semiisotropic
tau_p                    = 6.0  ;       parrinello-rahman is more stable with larger tau-p, DdJ, 20130422
compressibility          = 3e-4 3e-4
ref_p                    = 1.0 1.0
; Bonds and constraints
constraints              = none
constraint_algorithm     = Lincs
continuation             = yes
; Pull code
pull                    = yes
pull_ngroups            = 2
pull_ncoords            = 1
pull_group1_name        = LIP
pull_group2_name        = Protein
pull_coord1_type        = umbrella      ; harmonic biasing force
pull_coord1_groups      = 1 2
pull_coord1_rate        = 0.00002       ; in nm per ps so 0.001 is 1 nm per ns
pull_coord1_k           = 1000          ; kJ mol^-1 nm^-2
pull_coord1_start       = yes
pull-nstfout            = 300           ; write forces ca every 10 ps. nstxout-compressed should be the same
pull-nstxout            = 300           ; write com ca every 10 ps. nstxout-compressed should be the same
pull_coord1_geometry    = direction
pull-coord1-vec         = 0 0 1         ; z axis
EOF
gmx grompp -f pull.mdp -c ../eq.gro -r ../eq.gro -p ../topol_new.top -o pull.tpr -n ../index.ndx -quiet -maxwarn 1
gmx mdrun -v -deffnm pull -quiet -pin on -ntomp $noCPUs -ntmpi 1 -pinoffset $pinoffset

# for vizualization
echo 1 0 | gmx trjconv -f pull.xtc -s pull.tpr -pbc whole -center -dt 1000 -o pull_pbc.xtc -quiet

cat <<EOF > extract_frames.py
import numpy as np
# input file name
file_in = 'pull_pullx.xvg'
# output file name
file_out = 'extract_frames.ndx'
# generate output file
with open(file_out,'w') as f:
        f.write(' [ frames ] \n\n')
# find and add frames to output file
time,dist = np.genfromtxt(file_in,skip_header=17,skip_footer=3,unpack=True)
# skip the last few lines - might be corrupted as pull is stopped with error when dist reach half the box size
min_dist = dist[0]
max_dist = 7.0 # protein pulled far enough from membrane to be converged
step_size = 0.05 # stepsize in nm
dist_steps = np.arange(min_dist,max_dist,step_size)
index_prev = 0
for d in dist_steps:
    indices = np.where(dist>d)
    index = indices[0][0]
    frame = index + 1
    if index != index_prev:
        with open(file_out,'a') as f:
            f.write('%d\n' % frame)
    index_prev = index
with open(file_out,'a') as f:
        f.write('\n')
EOF
python extract_frames.py
echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -fr extract_frames.ndx -o extract_frames.xtc -quiet

## push
cp pull.mdp push.mdp
sed -i -e 's/pull_coord1_rate        = 0.00002/pull_coord1_rate        = -0.00002/g' push.mdp
gmx grompp -f push.mdp -c ../eq.gro -r ../eq.gro -p ../topol_new.top -o push.tpr -n ../index.ndx -quiet -maxwarn 1
gmx mdrun -v -deffnm push -quiet -pin on -ntomp $noCPUs -ntmpi 1 -pinoffset $pinoffset -nsteps 4000000 # 2.5 nm

# for vizualization
echo 1 0 | gmx trjconv -f push.xtc -s push.tpr -pbc whole -center -dt 1000 -o push_pbc.xtc -quiet

cp extract_frames.py extract_frames_push.py
sed -i -e "s/file_in = 'pull_pullx.xvg'/file_in = 'push_pullx.xvg'/g" extract_frames_push.py
sed -i -e "s/file_out = 'extract_frames.ndx'/file_out = 'extract_frames_push.ndx'/g" extract_frames_push.py
sed -i -e "s/min_dist = dist\[0\]/min_dist = dist[-1]/g" extract_frames_push.py
sed -i -e "s/max_dist = 7.0/max_dist = dist[0]/g" extract_frames_push.py
sed -i -e "s/np.where(dist>d)/np.where(dist<d)/g" extract_frames_push.py
python extract_frames_push.py
echo 0 | gmx trjconv -f push.xtc -s push.tpr -fr extract_frames_push.ndx -o extract_frames_push.xtc -quiet

cd ../..

fi

######################################################################################################
### GROMP AND SEND TO AXON ###########################################################################
######################################################################################################

if [ $GROMPP -eq 1 ]
then

cd ${protein}/umbrella

((Nframes_push=$(wc -l < extract_frames_push.ndx)-3))
((Nframes=$(wc -l < extract_frames.ndx)-3+$(wc -l < extract_frames_push.ndx)-3))
echo Nframes = $Nframes, Nframes_push = $Nframes_push

cat<<EOF > umbrella.mdp
integrator           = md
tinit                = 0.0
dt                   = 0.035                 ; 35 fs
nsteps               = 30000000              ; 30 * 35 ns = 1050 ns
nstxout              = 1000000
nstvout              = 1000000
nstfout              = 1000000
nstlog               = 10000
nstenergy            = 10000
nstxout-compressed   = 10000                 ; write every 350 ps
compressed-x-precision = 1000
nstlist              = 10
ns_type              = grid
pbc                  = xyz
;rlist                = 1.1
coulombtype          = Reaction_field
rcoulomb_switch      = 0.0
rcoulomb             = 1.1
epsilon_r            = 15
vdw_type             = cutoff
rvdw_switch          = 0.9
rvdw                 = 1.1
cutoff-scheme        = verlet
coulomb-modifier     = Potential-shift
vdw-modifier         = Potential-shift
epsilon_rf           = 0
verlet-buffer-tolerance  = 0.005
tcoupl               = v-rescale
tc-grps              = Protein LIP SOL
tau_t                = 1.0 1.0 1.0
ref_t                = 323 323 323
Pcoupl               = parrinello-rahman
Pcoupltype           = semiisotropic
tau_p                = 12.0
compressibility      = 3e-4 3e-4
ref_p                = 1.0 1.0
gen_vel              = no
gen_temp             = 323
gen_seed             = -1
constraints          = none
constraint_algorithm = Lincs
continuation         = no
lincs_order          = 4
lincs_warnangle      = 30
refcoord_scaling     = com
; pull code
pull                    = yes
pull_ngroups            = 2
pull_ncoords            = 1
pull_group1_name        = LIP
pull_group2_name        = Protein
pull_coord1_type        = umbrella      ; harmonic biasing force
pull_coord1_groups      = 1 2
pull_coord1_rate        = 0             ; don't pull, i.e. keep in position
pull_coord1_k           = 1000          ; kJ mol^-1 nm^-2
pull_coord1_start       = yes
pull-nstfout            = 50
pull_coord1_geometry    = direction
pull-coord1-vec         = 0 0 1         ; z axis
EOF

for frame in $(seq 1 $Nframes)
do 
  if [ $frame -le $Nframes_push ]
  then
    cat << EOF > frame.ndx
 [ frames ]

$frame

EOF
    echo 0 | gmx trjconv -f extract_frames_push.xtc -s push.tpr -fr frame.ndx -o frame.gro -quiet
  else
    cat << EOF > frame.ndx
 [ frames ]

$((frame-Nframes_push))

EOF
    echo 0 | gmx trjconv -f extract_frames.xtc -s pull.tpr -fr frame.ndx -o frame.gro -quiet
  fi
  echo -----------------------------------
  echo GROMPPing frame $frame of $Nframes
  echo -----------------------------------
  gmx grompp -f umbrella.mdp -c frame.gro -n ../index.ndx -p ../topol_new.top -o umbrella_step$frame.tpr -quiet -maxwarn 1
  rm frame.gro frame.ndx
done  

## send to axon 
mkdir -p $axon_dir
cp umbrella_step*.tpr $axon_dir

cat << EOF > Flow_axon.sh
#!/bin/bash 

## Define if you are in testing mode (on short queue) or running

## Set the job name.
#SBATCH --job-name=$protein

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
mv Flow_axon.sh $axon_dir
scp -r $axon_dir axon:/home/bioc1642/
rm -r $axon_dir

cd ../..

fi

######################################################################################################
### RUN LOCALLY ######################################################################################
######################################################################################################

if [ $RUN_LOCAL -eq 1 ]
then

cd ${protein}/umbrella

# run locally
#for frame in $(seq 1 $Nframes)
#do
#  gmx mdrun -v -deffnm umbrella_step$frame -quiet -pin on -ntopm $noCPUs -ntmpi 1 -pinoffset $pinoffset
#done

cd ../..

fi

######################################################################################################
### DOWNLOAD FROM AXON ###############################################################################
######################################################################################################

if [ $DOWNLOAD -eq 1 ]
then

cd ${protein}/umbrella

## download from axon
scp -r axon:/home/bioc1642/$axon_dir/*.tpr .
scp -r axon:/home/bioc1642/$axon_dir/*pullf.xvg .

cd ../..

fi

######################################################################################################
### ANALYZE WHAM #####################################################################################
######################################################################################################

if [ $ANALYZE -eq 1 ]
then

cd ${protein}/umbrella

echo analyze with wham

((Nframes=$(wc -l < extract_frames.ndx)-3+$(wc -l < extract_frames_push.ndx)-3))
echo $Nframes

#mkdir -p $axon_dir
#cd $axon_dir

#for frame in $(seq 1 $Nframes)
#do
#  echo "umbrella_step${frame}.tpr" >> tpr_files_umbrella.dat
#  echo "umbrella_step${frame}_pullf.xvg" >> pullf_files_umbrella.dat
#done

taskset -c 0-1  gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 323 -b 100000 -quiet -min 3 -max 7
taskset -c 0-1 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 323 -b 100000 -quiet -nBootstrap 200 -min 3 -max 7
#xmgrace bsProfs.xvg
#xmgrace bsResult.xvg

cd ../..

fi

######################################################################################################
### CLEAN UP #########################################################################################
######################################################################################################

## clean up
cd ${protein}/umbrella
rm -f \#*
cd ..
rm -f \#*
cd ..
rm -f \#*

