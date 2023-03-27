
###############################################
#################### PTEN #####################
###############################################

# go to folder for PTEN_FL selected binding mode
cd 4/PCPSPIP2_13/

## make gro file
echo 1 1 | gmx trjconv -s md_ref.tpr -f key_frame.gro -o key_frame_cent_prot.gro -pbc mol -center -ur compact -quiet

## remove C2 from gro file
cp key_frame_cent_prot.gro key_frame_cent_prot_Ptase.gro
# DO manually: remove all protein atoms from nubmer 187 to end. change total number of atoms to 444. Keep last line (box dimensions)

## change mdp file
cp /sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/production.mdp prod_Ptase.mdp
# DO manually: remove everything regarding LIP and SOLV (temperature coupling)

## change topology file
cp ../../partial_umbrella/Protein_PTEN_Pt.itp .# copy topology for PTEN_Ptase to folder
cp topol.top topol_Ptase.top
# DO manually: remove all lines from topology file molecule list, except Protein 1
# DO manually: change #include Protein.itp to #include Protein_PTEN_Pt.itp

## GROMPP gro file
gmx grompp -f prod_Ptase.mdp -c key_frame_cent_prot_Ptase.gro -r key_frame_cent_prot_Ptase.gro -p topol_Ptase.top -o key_frame_cent_prot_Ptase_from_PTEN_FL.tpr -quiet

## copy to PTEN_Ptase folder
cp key_frame_cent_prot_Ptase_from_PTEN_FL.tpr ../../3

## go to PTEN Ptase folder
cd ../../3

## RZZ analysis, with respect to Ptase from PTEN_FL binding site
for rep in {0..24}
do
    cd PCPSPIP2_$rep
    echo 1 | gmx rotmat -s ../key_frame_cent_prot_Ptase_from_PTEN_FL.tpr -f md_prot.xtc -fitxy -o Rzz.xvg -quiet
    cd ..
done

###############################################
#################### SHIP2 ####################
###############################################

# go to folder for SHIP2_FL selected binding mode
cd 5/PCPSPIP2_13/

## make gro file
echo 1 1 | gmx trjconv -s key_frame.tpr -f key_frame.gro -o key_frame_cent_prot.gro -pbc mol -center -ur compact -quiet

## remove C2 from gro file
cp key_frame_cent_prot.gro key_frame_cent_prot_Ptase.gro
# DO manually: remove all protein atoms before 3GlU and after 314SER. change total number of atoms to 728. Keep last line (box dimensions)

## change mdp file
cp /sansom/s157/bioc1642/Desktop/prj_C2/DLL4_MARTINI20/production.mdp prod_Ptase.mdp
# DO manually: remove everything regarding LIP and SOLV (temperature coupling)

## change topology file
cp ../../2/PCPSPIP2_0/Protein_B.itp Protein_SHIP2_Pt.itp # copy topology for SHIP2_Ptase to folder 
cp topol.top topol_Ptase.top
# DO manually: remove all lines from topology file molecule list, except Protein 1
# DO manually: change #include Protein.itp to #include Protein_SHIP2_Pt.itp
# DO manually: change Protein to Protein_B in molecule list

## GROMPP gro file
gmx grompp -f prod_Ptase.mdp -c key_frame_cent_prot_Ptase.gro -r key_frame_cent_prot_Ptase.gro -p topol_Ptase.top -o key_frame_cent_prot_Ptase_from_SHIP2_FL.tpr -quiet

## copy to SHIP2_Ptase folder
cp key_frame_cent_prot_Ptase_from_SHIP2_FL.tpr ../../2

## go to SHIP2_Ptase folder
cd ../../2

## RZZ analysis, with respect to Ptase from PTEN_FL binding site
for rep in {0..24}
do
    cd PCPSPIP2_$rep
    echo 1 | gmx rotmat -s ../key_frame_cent_prot_Ptase_from_SHIP2_FL.tpr -f md_prot.xtc -fitxy -o Rzz.xvg -quiet
    cd ..
done

