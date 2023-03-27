# load gmx etc. 
source ~/.bashrc

# user input
protein_folder=4 # 2=SHIP2_Ptase, 4=PTEN_FL2 (with N-terminal helix), 5=SHIP2_FL
lip=PCPSPIP2 # PC or PCPS or PCPSPIP2
rep=13 # SHIP2_Ptase (mode 2 from paper): rep 17, PTEN_FL2: rep=13, SHIP2_FL: rep=13 (no typo, it is rep=13 in both cases)

# performance
pinoffset=4

# go to CG2AT folder
cd /sansom/s157/bioc1642/Desktop/prj_PTEN_SHIP/${protein_folder}/${lip}_${rep}/CG2AT/FINAL

# pdb (output from cg2at)
if [ $protein_folder -eq 4 ]
then
    pdb=final_cg2at_de_novo.pdb
else
    pdb=final_cg2at_aligned.pdb
fi

#mdp files
mdp_folder=/sansom/s157/bioc1642/Desktop/prj_Notch/AT
min=$mdp_folder/min.mdp
nvt=$mdp_folder/nvt.mdp
npt=$mdp_folder/npt.mdp
prod=$mdp_folder/prod.mdp


# make index file
if [ $protein_folder -eq 4 ]
then
cat << EOF > input
r POPC POPS POPI
name 26 LIP
1 | 26
ri 184-352
name 28 C2
ri 1-183
name 29 Pt
q
EOF
elif [ $protein_folder -eq 2 ]
then
cat << EOF > input
r POPC POPS POPI
name 26 LIP
1 | 26
q
EOF
else
cat << EOF > input
r POPC POPS POPI
name 26 LIP
1 | 26
ri 324-461
name 28 C2
ri 1-323
name 29 Pt
q
EOF
fi
gmx make_ndx -f $pdb -quiet < input
rm input

# nvt equilibration
gmx grompp -f $nvt -c $pdb -r $pdb -p topol_final.top -n index.ndx -o nvt.tpr -quiet
gmx mdrun -deffnm nvt -v -quiet -pin on -pinoffset $pinoffset -ntomp 4 -ntmpi 1 -nb gpu

# calculate prot-lip COM distance of nvt
gmx distance -s nvt.tpr -f nvt.xtc -n index.ndx -tu ns -select "com of group LIP plus com of group Protein" -oall dist_nvt -quiet
gmx distance -s nvt.tpr -f nvt.xtc -n index.ndx -tu ns -select "com of group LIP plus com of group Protein" -oall dist_nvt -quiet

# npt equilibration
gmx grompp -f $npt -c nvt.gro -r nvt.gro -t nvt.cpt -p topol_final.top -n index.ndx -o npt.tpr -quiet
gmx mdrun -deffnm npt -v -quiet -pin on -pinoffset $pinoffset -ntomp 4 -ntmpi 1 -nb gpu

# prepare and test production run
gmx grompp -f $prod -c npt.gro -r npt.gro -t npt.cpt -p topol_final.top -n index.ndx -o prod.tpr -quiet
gmx mdrun -deffnm prod -v -quiet -pin on -pinoffset $pinoffset -ntomp 4 -ntmpi 1 -nb gpu -nsteps 10000

# production run
#gmx mdrun -deffnm prod -v -quiet -pin on -pinoffset 0 -ntomp 4 -ntmpi 1 -nb auto -nsteps 20000000 # 40 ns
#gmx mdrun -deffnm prod -v -quiet -pin on -pinoffset 0 -ntomp 4 -ntmpi 1 -nb auto -nsteps 20000000 # 40 ns

# extend by 200,000 ps = 200 ns (and append to previous files)
#gmx convert-tpr -s prod.tpr -extend 200000 -o prod.tpr -quiet

# continue production run
cat << EOF > continue.sh
gmx mdrun -deffnm prod -cpi prod.cpt -append -v -quiet -pin on -pinoffset $pinoffset -ntomp 4 -ntmpi 1 -nb gpu
EOF
nohup bash continue.sh > nohup_folder_${protein_folder}_AT.out &

# center protein in box, make sparser trajectory (from 10 ps (600 ns -> 60000 frames) to 1000 ps per frame (600 ns -> 600 frames))
echo 1 0 | gmx trjconv -s prod.tpr -f prod.xtc -o prod_cent.xtc -pbc mol -center -ur compact -dt 1000 -quiet

## Analyze result

# calculate rmsd with respect to CG frame
echo 1 1 | gmx rms -s npt.tpr -f prod_cent.xtc -n index.ndx -tu ns -o rmsd.xvg -quiet

# calculate prot-lip COM distance of prod
gmx distance -s prod.tpr -f prod_cent.xtc -n index.ndx -tu ns -select "com of group LIP plus com of group Protein" -oall dist -quiet
gmx distance -s prod.tpr -f prod_cent.xtc -n index.ndx -tu ns -select "com of group LIP plus com of group C2" -oall dist_C2 -quiet
gmx distance -s prod.tpr -f prod_cent.xtc -n index.ndx -tu ns -select "com of group LIP plus com of group Pt" -oall dist_Pt -quiet

# make index group for C2
gmx grompp -f $nvt -c $pdb -r $pdb -p topol_final.top -n index.ndx -o nvt.tpr -quiet

# center first and last frame, for vizualization
echo 1 27 | gmx trjconv -s npt.tpr -f npt.gro -n index.ndx -o first.gro -center -pbc mol -quiet
echo 1 27 | gmx trjconv -s prod.tpr -f prod.gro -n index.ndx -o last_cent.gro -center -pbc mol -quiet
echo 27 | gmx trjconv -s prod.tpr -f prod.gro -n index.ndx -o last.gro -pbc mol -quiet
echo 1 0 | gmx trjconv -s first.gro -f last.gro -o last_fit.gro -fit rotxy+transxy -quiet

# calculate C2-Pt COM distance of prod
gmx distance -s prod.tpr -f prod_cent.xtc -n index.ndx -tu ns -select "com of group C2 plus com of group Pt" -oall dist_C2-Pt -quiet

# calculate prot-PO2 COM distance of prod
if [ $protein_folder -eq 4 ]
then
imin=778
imax=800
elif [ $protein_folder -eq 5 ]
then
imin=888
imax=909	
fi

for i in $(seq $imin $imax) 
do
echo i = $i
cat << EOF > input
ri $i
name 30 PI
q
EOF
gmx make_ndx -f $pdb -n index -o index2.ndx -quiet < input
rm input
echo 1 30 | gmx mindist -s prod.tpr -f prod_cent.xtc -n index2.ndx -tu ns -od mindist_Prot-PI$i -quiet
echo 28 30 | gmx mindist -s prod.tpr -f prod_cent.xtc -n index2.ndx -tu ns -od mindist_C2-PI$i -quiet
echo 29 30 | gmx mindist -s prod.tpr -f prod_cent.xtc -n index2.ndx -tu ns -od mindist_Pt-PI$i -quiet
rm index2.ndx
done

# clean up
rm -f \#*

# back to parent directory
cd ../../../../
