#!/usr/bin/env python
import sys,os
import matplotlib.pyplot as plt
from wavesolver_results import *
from ConfigParser import SafeConfigParser,NoOptionError
import scipy
from scipy import signal

if len(sys.argv) != 3:
    print '**Usage: %s <data_folder> <audio_filename>' %(sys.argv[0])
    sys.exit()

results = Wavesolver_Results()
results.Set_Folder(sys.argv[1])
data = results.Read_All_Audio(False, 1, None, None, filename=sys.argv[2])
print data.shape

plt.figure()
plt.plot(data)
plt.show()
