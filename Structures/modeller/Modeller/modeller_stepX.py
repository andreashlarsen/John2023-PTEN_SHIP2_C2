from modeller import *
from modeller.automodel import *

log.verbose()    # request verbose output

class MyModel(automodel):
    def select_atoms(self):
        return selection(self) - selection(self.residue_range('22:', '280:'),self.residue_range('314:', '350:'))

    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
        rsr.add(secondary_structure.alpha(self.residue_range('1:', '11:')))

env = environ()
env.io.atom_files_directory = ['.', '../atom_files']

a = MyModel(env, alnfile='work-2.ali', knowns=('1D5R_modified'), sequence='PTEN-seq')
a.starting_model = 1
a.ending_model = 1
a.make()
