#!/usr/bin/env python

import my_io as IO
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm 
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

# writing the wav files
for ii in range(N_points): 
    scipy.io.wavfile.write('point_%u.wav' %(ii), stepRate, listenedData[:,ii]/maxValue)


listenedData = listenedData[:, ::10]
listeningPositions = listeningPositions[::10, :]
N_points = listeningPositions.shape[0]
print listenedData.shape, listeningPositions.shape

plotting = True
if plotting:
    # plotting
    colors = [cm.jet(x) for x in np.linspace(1, 0, N_points)]
    fig = plt.figure(figsize=[10, 10])
    for ii in range(N_points): 
        if ii == 0 or ii == N_points-1:
            label = '%.2f %.2f %.2f' %(listeningPositions[ii, 0], listeningPositions[ii, 1], listeningPositions[ii, 2])
            plt.plot(listenedData[:,ii], label=label, color=colors[ii], linewidth=2) 
        else:
            plt.plot(listenedData[:,ii], color=colors[ii], linewidth=2) 
    plt.legend(loc=4)
    plt.xlabel('frame') 
    plt.ylabel('Pressure (Pascal)')
    plt.title('Line samples pressure time series (RAS)')
    plt.xlim([0, 500])
    plt.grid()
    # plt.title('Monitor points for ball drop using ghost cell')
    
    if len(sys.argv) == 5: 
        plot_index = int(sys.argv[4]) 
        data = listenedData[:, plot_index] 
        plt.figure() 
        plt.plot(data, label=listeningPositions[plot_index, :]) 
        plt.xlabel('frame') 
        plt.ylabel('Pressure (Pascal)')
        plt.xlim([0, 400])
        plt.grid()
        print 'Listening point = ', listeningPositions[plot_index, :]
        # try to extract the max in the last couple of cycles
        extractEndPercentage = 0.2
        N_extract = N_steps * extractEndPercentage
        dataExtracted = data[-N_extract::]
        print 'End %.2f percent of the listened data was extracted. Its max = %f' %(extractEndPercentage*100.0, np.absolute(dataExtracted).max())
    
    plt.show()
