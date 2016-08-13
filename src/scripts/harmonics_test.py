#!/usr/bin/env python
import sys
import numpy as np 
import matplotlib.pyplot as plt

################################################################################
################################################################################
def CompareVelocityBC():
    fileFBem = 'allVelocityFBem.txt'
    fileFDTD = 'allVelocityFDTD.txt'

    dataFBem = np.loadtxt(fileFBem)
    dataFDTD = np.loadtxt(fileFDTD)

    plt.figure() 
    plt.plot(dataFBem[:,1], label='FBem')
    plt.plot(dataFDTD[:,1], label='FDTD')
    plt.legend()
    plt.show()

################################################################################
################################################################################
def PlotHarmonicsOutput():
    if len(sys.argv) != 2: 
        print '**Usage: %s <file>' %(sys.argv[0])
        sys.exit()
    
    filename = sys.argv[1]
    print filename
    data = np.loadtxt(filename)
    print data
    
    
    plt.figure() 
    plt.plot(data[:, 0], data[:, 1], '-o')
    plt.show()

################################################################################
################################################################################
if __name__ == '__main__':
    # CompareVelocityBC()
    PlotHarmonicsOutput()
