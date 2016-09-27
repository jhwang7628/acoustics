#!/usr/bin/env python 
import my_io as io
import scipy
import matplotlib.pyplot as plt 
import sys 

if len(sys.argv) != 3: 
    print '**Usage: %s <ODE_result_displacement_file> <mode>' %(sys.argv[0])
    sys.exit()

mode = int(sys.argv[2])
if mode == 0: 
    data = io.readMatrixXdBinary(sys.argv[1])
    
    # samples = range(0, data.shape[1], 1000)
    samples = [49757]
    for sampleVertexID in samples: 
        # sampleVertexID = 49757
        vertexDisplacement = data[:,sampleVertexID]
        wavname = sys.argv[1] + '.%u.wav' %(sampleVertexID)
        scipy.io.wavfile.write(wavname, 40000, vertexDisplacement / max(abs(vertexDisplacement)))
    
    plt.figure()
    plt.plot(vertexDisplacement, '-o') 
    plt.show()

elif mode == 1: 
    data = io.readMatrixXdBinary(sys.argv[1])
    print data
    
    samples = range(0, data.shape[1])
    for sampleVertexID in samples: 
        # sampleVertexID = 49757
        vertexDisplacement = data[:,sampleVertexID]
    
    plt.figure()
    plt.plot(vertexDisplacement, '-o') 
    plt.show()
