# 
# MODELLER LOG
#
# use scripts from: 
# https://salilab.org/modeller/wiki/Missing%20residues
#
# IMPORT MODELLER
module load modeller/9v16/64
#
# STEP 0
#name=5okm_FL
#cp /sansom/s157/bioc1642/Desktop/prj_C2/Structures/${name}.pdb .
cp ../../1d5r.pdb .
sed -i '/HETATM/d' 1d5r.pdb
mv 1d5r.pdb 1D5R_modified.pdb
#
# STEP 1
#
# modify modeller_step1.py
#mod9.16 modeller_step1.py
# this script generates *.seq
#
# now generate alignment file 
#cp *.seq alignment.ali
# follow instructions on webpage (and use missing residue info in pdb file) to make alignment.ali
# https://salilab.org/modeller/wiki/Missing%20residues
#
# STEP 2
# modify modeller_step2.py
mod9.16 modeller_stepX.py
#
# Choose best model
#cp ${name}_fill.B99990001.pdb /sansom/s157/bioc1642/Desktop/prj_C2/Structures/5okm_FL_clean.pdb
