#!/usr/bin/env python

import sys
import numpy as np
from scipy.io import wavfile
import glob
import matplotlib.pyplot as plt

files = sorted(glob.glob('*.wav'))
print 'read files: ', files
data = []

for f in files: 
    print f
    r, d = wavfile.read(f)
    d = d / float(2**15)
    localData = dict()
    localData['name'] = f
    localData['rate'] = r
    localData['signal'] = d
    data.append(localData)

N_plot = 7
plt.figure()
for d in range(N_plot):
    filename = data[d]['name']
    if filename.find('synthesized') != -1: 
        if filename.find('attenuated') == -1:
            color = 'r'
        else: 
            # color = '#ff7f7f'
            color = '#b20000'
    else: 
        color = 'b'

    N_data = len(data[d]['signal'])
    dt = 1./float(data[d]['rate'])
    t = np.linspace(0, float(N_data)*dt, N_data)
    plt.subplot(N_plot, 1, d+1)
    plt.plot(t, data[d]['signal'], color=color)
    plt.ylim([-1, 1])
    plt.grid()
    plt.title(filename)

    if (d != N_plot-1): 
        plt.tick_params(labelbottom='off')

plt.xlabel('time (s)')
plt.show()
