#!/usr/bin/env python
import sys,os
import matplotlib.pyplot as plt
from wavesolver_results import *
from ConfigParser import SafeConfigParser,NoOptionError
import scipy
from scipy import signal

if len(sys.argv) < 3:
    print '**Usage: %s <data_folder> <audio_filename> [frame_rate] [offset_time]' %(sys.argv[0])
    sys.exit()

results = Wavesolver_Results()
results.Set_Folder(sys.argv[1])
data = results.Read_All_Audio(False, 1, None, None, filename=sys.argv[2], N_speakers_along_ray=30)

if len(sys.argv) > 3:
    dt = 1./float(sys.argv[3])
    if len(sys.argv) > 4:
        t0 = float(sys.argv[4])
    else:
        t0 = 0.0
    t = np.linspace(t0, t0+dt*data.shape[0], data.shape[0])
else:
    t = range(data.shape[0])

plt.figure()
plt.plot(t, data)
plt.show()

sampfreq = 1./dt
wavfreq  = 44100
prefix   = 'point'
wavformat= '16int'
out_dir  = 'post'
N_points = data.shape[1]
rateRatio = float(sampfreq) / float(wavfreq)
for ii in range(N_points):
    print 'point %u' %(ii)
    # scipy resample uses FFT, zero-pad it to make it fast
    outputdata = data[:,ii].copy()
    nextpow2 = int(np.power(2, np.floor(np.log2(len(outputdata)))+1))
    T = float(outputdata.shape[0])/float(sampfreq)
    outputdata = np.pad(outputdata, (0, nextpow2-len(outputdata)), 'constant', constant_values=(0.,0.))
    outputdata = signal.resample(outputdata, int(float(len(outputdata))/rateRatio))
    outputdata = outputdata[:int(round(T*wavfreq))]
    normalization = np.absolute(outputdata).max()
    if wavformat == '32float':
        finaldata = outputdata/normalization
    elif wavformat == '16int':
        finaldata = ((outputdata/normalization*32767)).astype('int16')
    scipy.io.wavfile.write('%s/%s_%u.wav' %(out_dir, prefix, ii), wavfreq, finaldata)

