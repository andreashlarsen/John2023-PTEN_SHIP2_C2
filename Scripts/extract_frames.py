#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This program look in the output file from steered md
1) finds the frames with given distances (protein-membrane)
2) writes frame numbers into an index file that can be used by gmx trjconv

"""

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
max_dist = 7.0       # protein pulled far enough from membrane to be converged
step_size = 5     # in hundreds of nm (to get an integer and avoid floats in bash script), use this step size of to get good histogram overlap for WHAM
step_size_nm = step_size/100.0 # stepsize in nm
dist_steps = np.arange(min_dist,max_dist,step_size_nm)

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
