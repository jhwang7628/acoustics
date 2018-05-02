#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from signal_processing import *
from wavesolver_results import *
from scipy.interpolate import interp1d

################################################################################
################################################################################
def ReadModes(txt, mode_ids):
    mode_ids_set = set(mode_ids)
    freq_map = None
    with open(txt,'r') as stream:
        freq_map = dict()
        lines = stream.readlines()
        for l in lines:
            tokens = l.split()
            if len(tokens) >= 6 and int(tokens[1]) in mode_ids_set:
                freq_map[int(tokens[1])] = float(tokens[4])
    return freqa_map

################################################################################
################################################################################
data_dir = sys.argv[1]
pt = 1
N_window = 1024
dN = 512
samp_freq = 118826.
dt = 1./samp_freq
# freqs = [489.4, 2192.8, 3076.7, 3373., 4013.9, 4447.5, 7152.0]
freqs = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390]
skip = 5

results = Wavesolver_Results()
results.Set_Folder(data_dir)
all_data = results.Read_All_Audio(False, 1, 1, N_speakers_along_ray=6)
N = all_data.shape[0]
YF = np.zeros((N/dN,N_window/2))
fftdata = np.zeros((N/dN, len(freqs)))

plt.figure()
for f in freqs:
    plt.loglog([f,f],[1,1E11],'--k')
colors = [cm.jet(x) for x in np.linspace(0, 1, N/dN)]
for ii in range(N/dN):
    print 'computing %u/%u' %(ii, N/dN)
    if (ii*dN+N_window) > N:
        break
    xf,yf,ymag = ComputeAverageFFT(dt, all_data[ii*dN:ii*dN+N_window,pt], N=N_window)
    if (ii % skip == 0):
        plt.loglog(xf, ymag, color=colors[ii])

    YF[ii,:] = ymag
    f2 = interp1d(xf, ymag, kind='cubic')
    ymag_interp = f2(freqs)
    fftdata[ii,:] = ymag_interp

plt.figure()
colors = [cm.jet(x) for x in np.linspace(0, 1, len(freqs))]
for ii,f in enumerate(freqs):
    norm = sum(abs(fftdata[10:-1,ii]))*float(len(fftdata[10:-1,ii]))
    plt.plot(fftdata[:-1,ii]/norm, '-o', color=colors[ii], label=f)
plt.legend(loc=4)
plt.figure()
plt.plot(YF)
plt.show()
