#!/usr/bin/env python 
import sys
import numpy as np 
import matplotlib.pyplot as plt

################################################################################
################################################################################
def ReadCommentLine(inFile): 
    buf = inFile.readline()
    if buf.find('# ') == -1: 
        raise IOError('Not comment line: %s' %(buf))

################################################################################
################################################################################
if __name__ == '__main__': 
    if len(sys.argv) != 2: 
        print '**Usage: %s <filename>' %(sys.argv[0])
        sys.exit()
    
    # read file in buffer
    filename = sys.argv[1]
    inFile = open(filename)
    
    # start parse
    ReadCommentLine(inFile)
    buf = inFile.readline().split()
    N_slices = int(buf[0])
    sliceDataPointer = int(buf[1])
    slices = dict()
    for s_idx in range(N_slices): 
        ReadCommentLine(inFile)
        buf = inFile.readline().split()
        sliceIndex = int(buf[0]) 
        N_vertices = int(buf[1]) 
        N_dataDim  = int(buf[2]) 
        ReadCommentLine(inFile)
        data = np.zeros((N_vertices, N_dataDim+3))
        for v_idx in range(N_vertices): 
            buf = inFile.readline().split()
            data[v_idx, 0] = float(buf[0]) 
            data[v_idx, 1] = float(buf[1]) 
            data[v_idx, 2] = float(buf[2]) 
            for d_idx in range(N_dataDim): 
                data[v_idx, 3+d_idx] = float(buf[3+d_idx])
        slices[s_idx] = data

    # post-process
    plt.figure() 
    for s_idx in slices: 
        plt.plot(slices[s_idx][:, 3], '-o', label='normalized residual')
    plt.show()
