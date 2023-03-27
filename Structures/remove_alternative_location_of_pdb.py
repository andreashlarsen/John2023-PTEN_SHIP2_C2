import numpy as np

f = open('5okm_B_ptase.pdb','r')
ff = open('5okm_B_ptase_clean.pdb','w')

line = f.readline()
while line:
    if len(line) < 16:
        ff.write('%s' % line)
    elif line[16] == 'A':
        ff.write('%s %s' % (line[0:16],line[17:]))
    elif line[16] == 'B':
        pass
    else:
        ff.write('%s' % line)
    line = f.readline()

f.close()
ff.close()
