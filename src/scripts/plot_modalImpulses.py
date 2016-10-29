#!/usr/bin/env python 
import sys,math
import matplotlib.pyplot as plt
if len(sys.argv) != 2: 
    print '**Usage: %s <modal_impulse_file>' %(sys.argv[0]) 
    sys.exit(1)

filename = sys.argv[1]
lines = open(filename, 'r').readlines()
impactSpeed = []
time = []
for l in lines: 
    tokens = l.split()
    time.append(float(tokens[0]))
    impactSpeed.append(abs(float(tokens[3])))

plt.figure() 
plt.plot(time, impactSpeed) 
plt.show()
