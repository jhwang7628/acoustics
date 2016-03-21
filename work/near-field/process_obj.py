#!/usr/bin/env python 
import sys
import numpy as np

if len(sys.argv) < 3: 
    print '**Usage: %s <mesh_obj> <out_obj> [scale]' 
    sys.exit()


filename = sys.argv[1]
# outname = filename[-4]+'_modified.obj'
outname = sys.argv[2]

if len(sys.argv) == 4: 
    scale = float(sys.argv[3])
else: 
    scale = 1.0

lines = open(filename,'r').readlines()
centerPosition = np.zeros(3)
N_vertex = 0
obj = []
print 'scanning through the files to collect vertex position'
for l in lines:
    line = l.split(' ')
    if line[0] == 'v': 
        centerPosition[0] += float(line[1])
        centerPosition[1] += float(line[2])
        centerPosition[2] += float(line[3])
        N_vertex += 1

centerPosition /= float(N_vertex) 
print ' computed center=', centerPosition
print 'recentering'


for l in lines: 
    line = l.split(' ') 
    if line[0] == 'v': 
        obj.append('v %.16f %.16f %.16f\n' %(scale*float(line[1]) - centerPosition[0], scale*float(line[2]) - centerPosition[1], scale*float(line[3]) - centerPosition[2]))
    elif line[0] == 'f': 
        obj.append(l)


of = open(outname, 'w')
for l in obj: 
    of.write(l)
of.close()
