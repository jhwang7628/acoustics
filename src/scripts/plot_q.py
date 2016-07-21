#!/usr/bin/env python 
import my_io as io
import scipy
import matplotlib.pyplot as plt 
import sys 

if len(sys.argv) != 3: 
    print '**Usage: %s <ODE_result_q_file> <mode_index: -1 for all modes>' %(sys.argv[0])
    sys.exit()

data = io.readMatrixXdBinary(sys.argv[1])
mode = int(sys.argv[2])
if mode >= 0 and mode < data.shape[1]: 
    plt.figure()
    plt.plot(data[:,mode], '-o') 
    plt.show()
elif mode == -1:  # all modes
    plt.figure()
    plt.plot(data[:,:], '-o') 
    plt.show()
else: 
    print '%u mode out of bounds (0-%u)' %(mode, data.shape[1]-1)
