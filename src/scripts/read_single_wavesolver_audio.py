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
data = results.Read_All_Audio(False, 1, None, None, filename=sys.argv[2])
print data.shape

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
