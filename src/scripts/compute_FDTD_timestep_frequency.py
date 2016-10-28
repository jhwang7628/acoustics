#!/usr/bin/env python 
import sys,math

if len(sys.argv) < 2: 
    print '**Usage: %s <cell_size> [dimension=3]' %(sys.argv[0])
    sys.exit(1)

h = float(sys.argv[1])
if len(sys.argv) == 3: 
    d = int(sys.argv[2])
else: 
    d = 3

tmax = h / (math.sqrt(float(d))*343.)
print 'Min sampling frequency = %u (max time step size = %.12f), if h = %f for %u dimension' %(int(1./tmax), tmax, h, d)
