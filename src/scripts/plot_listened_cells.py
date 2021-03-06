#!/usr/bin/env python
import binary_io as IO
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm 
import sys
import glob
import scipy
from scipy import signal

def PadZero2D(data, N_front_pad, N_end_pad):
    N_rows = data.shape[0]
    N_cols = data.shape[1]
    newData = np.zeros((N_front_pad+N_rows+N_end_pad, N_cols))
    newData[N_front_pad:N_front_pad+N_rows, :] = data
    # newData[:, N_front_pad:N_front_pad+N_cols] = data

    return newData

def PadZero(data, N_front_pad, N_end_pad):
    N_rows = len(data)
    newData = np.zeros(N_front_pad+N_rows+N_end_pad)
    newData[N_front_pad:N_front_pad+N_rows] = data
    # newData[:, N_front_pad:N_front_pad+N_cols] = data

    return newData

if __name__ == '__main__':
    if (len(sys.argv) < 4): 
        print '**Usage: %s <listening_position_dat> <listening_data_dat_prefix> <rate(Hz)> [single_plot_cell]' %(sys.argv[0])
        sys.exit()
    
    plt.rcParams['agg.path.chunksize'] = 20000
    
    listeningPositions = IO.readMatrixXdBinary(sys.argv[1]) 
    filenames = sorted(glob.glob(sys.argv[2]+'[0-9]*'))
    N_steps = len(filenames)
    N_points = listeningPositions.shape[0]
    stepRate = int(sys.argv[3])
    wavRate = 44100.
    rateRatio = float(stepRate) / float(wavRate)
    plotting = True
    logPlot = False
    
    verifyAnalytical = False
    if verifyAnalytical: 
        print 'Verify with analytical solution'
        analyticalOutputFolder = 'analytical'
        analyticalFilenames = sorted(glob.glob(analyticalOutputFolder+'/'+sys.argv[2]+'[0-9]*'))
        analyticalListenedData = np.zeros((N_steps, N_points))
    print 'Loading %u files: ' %(len(filenames)), filenames[0], ',', filenames[1], ', ... ,' , filenames[-1]
    
    # get listened data for each step
    listenedData = np.zeros((N_steps, N_points))
    step = 0
    for f in filenames: 
        stepData = IO.readMatrixXdBinary(f, verbose=0)
        listenedData[step, :] = stepData.transpose()
        if (verifyAnalytical):
            stepData = IO.readMatrixXdBinary(analyticalFilenames[step], verbose=0)
            analyticalListenedData[step, :] = stepData.transpose()
    
        step += 1
    
    maxValue = np.absolute(listenedData).max()
    print 'Normalize all data by max value = %f' %(maxValue)
    # writing the wav files
    frontPadTime = 0.0
    endPadTime = 0.0
    N_frontPad = int(frontPadTime*float(wavRate))
    N_endPad = int(endPadTime*float(wavRate))
    N_frontPadStep = int(frontPadTime*float(stepRate))
    N_endPadStep = int(endPadTime*float(stepRate))
    allSound=np.array([0.0])
    # N_frontPad = 0
    for ii in range(N_points): 
        if len(sys.argv) == 5 and ii != int(sys.argv[4]): 
            continue
        print 'Creating wav file for listening position', listeningPositions[ii, :]
        outputData = listenedData[:, ii]
        # normalizationConstant = maxValue
        normalizationConstant = np.absolute(outputData).max()
        # if normalizationConstant > 1E-14:
        #     outputData /= normalizationConstant
        outputData = signal.resample(outputData, int(float(len(outputData))/rateRatio))
        outputData = PadZero(outputData.copy(), N_frontPad, N_endPad)
        #if (ii < 25):
        #    allSound = np.concatenate([allSound, outputData])
        scipy.io.wavfile.write('point_%u.wav' %(ii), wavRate, outputData/normalizationConstant)
    #scipy.io.wavfile.write('allpoints_%u.wav' %(ii), wavRate, allSound)
    
    # extract minimum and compare with 1/r decay 
    dataMin = np.zeros(N_points)
    dataR   = np.zeros(N_points)
    for pt_idx in range(N_points): 
        dataMin[pt_idx] = abs(min(listenedData[:, pt_idx]))
        dataR[pt_idx] = np.linalg.norm(listeningPositions[pt_idx, :])
    one_over_r = np.divide(np.ones(N_points), dataR)
   
    IO.writeMatrixXdBinary('test_listened_data.dat', listenedData)
    # trim data 
    # listenedData = listenedData[:, ::20]
    # listeningPositions = listeningPositions[::20, :]
    # N_points = listeningPositions.shape[0]
    
    if plotting:
        t = np.linspace(0., float(N_steps)*(1./float(stepRate)), N_steps)
        # plotting
        colors = [cm.jet(x) for x in np.linspace(1, 0, N_points)]
        lw = 2.0
        ls = '-'
        fig = plt.figure(figsize=[10, 10])
        if len(sys.argv) == 5: 
            plot_index = int(sys.argv[4]) 
            data = listenedData[:, plot_index]
            data = PadZero(data.copy(), N_frontPadStep, N_endPadStep)
            tPad = np.linspace(0., float(len(data))*(1./float(stepRate)), len(data))
            plt.plot(tPad, data, ls, label=listeningPositions[plot_index, :]) 
            #plt.plot(data, '-', label=listeningPositions[plot_index, :]) 
            plt.grid()
            if (verifyAnalytical):
                plt.plot(analyticalListenedData[:, plot_index])
                  
            print 'Listening point = ', listeningPositions[plot_index, :]
            # try to extract the max in the last couple of cycles
            extractEndPercentage = 0.2
            N_extract = N_steps * extractEndPercentage
            dataExtracted = data[-N_extract::]
            print 'End %.2f percent of the listened data was extracted. Its max = %f' %(extractEndPercentage*100.0, np.absolute(dataExtracted).max())
        else:
            for ii in range(N_points): 
                # if ii == 0 or ii == N_points-1:
                label = '%.4f %.4f %.4f' %(listeningPositions[ii, 0], listeningPositions[ii, 1], listeningPositions[ii, 2])
                if not logPlot:
                    plt.plot(listenedData[:,ii], label=label, color=colors[ii], linewidth=lw) 
                else: 
                    plt.semilogy(np.multiply(listenedData[:,ii], listenedData[:,ii]), label=label, color=colors[ii], linewidth=lw) 
                if verifyAnalytical:
                    print 'max at %u is %f' %(ii, max(analyticalListenedData[:, ii]))
                    plt.plot(analyticalListenedData[:, ii], 'x', color=colors[ii], )
                # else:
                    # plt.plot(listenedData[:,ii], color=colors[ii], linewidth=lw) 
            # plt.legend(loc=4)
        # plt.title('Collocated formulation with PML')
        plt.legend()
        plt.xlabel('frame') 
        plt.ylabel('Pressure (Pascal)')
        plt.grid()
        #plt.ylim([-2E-7, 2E-7])
        # plt.title('Monitor points for ball drop using ghost cell')
        
    
        # # plot 1/r 
        # fig = plt.figure(figsize=[10, 10])
        # plt.loglog(dataR, dataMin, label='dipole solution', linewidth=lw)
        # plt.loglog(dataR, one_over_r*30., 'k--', label='1/r decay rate', linewidth=lw)
        # plt.legend()
        # yticks = np.logspace(1, 3, 10)
        # xticks = np.logspace(-2, 0, 3)
        # plt.yticks(yticks, ['%.2e' %(x) for x in yticks] )
        # plt.xticks(xticks)
        # plt.xlabel('r (m)')
        # plt.ylabel('pressure (Pa)')
        # plt.title('1/r decay')
        # plt.axis('tight')
        
        plt.show()
