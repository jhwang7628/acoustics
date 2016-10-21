#!/usr/bin/env python
import numpy as np 
import sys

if len(sys.argv) != 3: 
    print '**Usage %s <old_objfile> <new_objfile>' %(sys.argv[0])
    sys.exit(1)

old_filename = sys.argv[1]
new_filename = sys.argv[2]
inFile = open(old_filename)
lines = inFile.readlines()
s = """"""
center = np.zeros(3) 
N_v = 0
for l in lines:
    if l.split(' ')[0] == 'v': 
        for dim in range(3): 
            center[dim] += float(l.split(' ')[1:4][dim])
        N_v += 1

center /= float(N_v)
for l in lines: 
    if l.split(' ')[0] == 'v': 
        newPoint = np.zeros(3)
        for dim in range(3): 
            newPoint[dim] = float(l.split(' ')[1:4][dim]) - center[dim]
        s += 'v %.16f %.16f %.16f\n' %(newPoint[0], newPoint[1], newPoint[2])
    else: 
        s += l

open(new_filename, 'w').write(s)

