#!/usr/bin/env python

import my_io as IO
import numpy as np 
import matplotlib.pyplot as plt
import sys
import glob
import scipy

if (len(sys.argv) != 3): 
    print '**Usage: %s <listening_position_dat> <listening_data_dat_prefix>' %(sys.argv[0])
    sys.exit()

listeningPositions = IO.readMatrixXdBinary(sys.argv[1]) 
filenames = sorted(glob.glob(sys.argv[2]+'*'))
N_steps = len(filenames)
N_points = listeningPositions.shape[0]

print 'Listening position: \n', listeningPositions
print 'Loaded files: ', filenames[0], ',', filenames[1], ', ... ,' , filenames[-1]

# get listened data for each step
listenedData = np.zeros((N_steps, N_points))
step = 0
for f in filenames: 
    stepData = IO.readMatrixXdBinary(f, verbose=0)
    listenedData[step, :] = stepData.transpose()
    step += 1

# plotting
plt.figure()
for ii in range(N_points): 
    plt.plot(listenedData[:,ii], label=listeningPositions[ii, :]) 
    scipy.io.wavfile.write('point_%u.wav' %(ii), 176400, listenedData[:,ii])
plt.legend()
plt.show()

# write wav file

