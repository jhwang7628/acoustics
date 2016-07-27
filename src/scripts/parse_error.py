#!/usr/bin/env python 
import sys 
import numpy as np 

lines = open(sys.argv[1]).readlines()
# points = np.zeros((len(lines),3))
points = []
outFile = open(sys.argv[2], 'w')

count=0
for l in lines: 
    a = l.find('[')
    b = l.find(']')
    point_s = l[a+1:b]
    outFile.write(point_s)
    outFile.write('\n')
    # points.append(point_s)

# print points

