#!/usr/bin/env python

import my_io as IO
import numpy as np 
import matplotlib.pyplot as plt
import sys
import glob
import scipy

if (len(sys.argv) < 4): 
    print '**Usage: %s <listening_position_dat> <listening_data_dat_prefix> <rate(Hz)> [single_plot_cell]' %(sys.argv[0])
    sys.exit()

listeningPositions = IO.readMatrixXdBinary(sys.argv[1]) 
filenames = sorted(glob.glob(sys.argv[2]+'[0-9]*'))
N_steps = len(filenames)
N_points = listeningPositions.shape[0]
stepRate = int(sys.argv[3])

print 'Listening position: \n', listeningPositions
print 'Loading %u files: ' %(len(filenames)), filenames[0], ',', filenames[1], ', ... ,' , filenames[-1]

# get listened data for each step
listenedData = np.zeros((N_steps, N_points))
step = 0
for f in filenames: 
    stepData = IO.readMatrixXdBinary(f, verbose=0)
    listenedData[step, :] = stepData.transpose()
    step += 1

maxValue = np.absolute(listenedData).max()
print 'Normalize all data by max value = %f' %(maxValue)

# plotting
plt.figure()
for ii in range(N_points): 
    plt.plot(listenedData[:,ii]/maxValue, label=listeningPositions[ii, :]) 
    scipy.io.wavfile.write('point_%u.wav' %(ii), stepRate, listenedData[:,ii]/maxValue)
plt.legend(loc=2)
plt.xlabel('frame') 
plt.ylabel('Pressure (Pascal)')
# plt.title('Monitor points for ball drop using ghost cell')

if len(sys.argv) == 5: 
    plot_index = int(sys.argv[4]) 
    data = listenedData[:, plot_index] 
    plt.figure() 
    plt.plot(data[::4], label=listeningPositions[plot_index, :]) 
    plt.xlabel('frame') 
    plt.ylabel('Pressure (Pascal)')

plt.show()


