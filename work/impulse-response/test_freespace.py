#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt 
import sys 
from SignalProcessing import ComputeSpectrum


if len(sys.argv) != 2: 
    print '**Usage: %s <interpolated_data>' %(sys.argv[0])
    sys.exit()

l = np.loadtxt('head_listening.txt') 

a = np.loadtxt(sys.argv[1]) 

colormap = plt.cm.jet
num_plots_1 = a.shape[1]
# num_plots_1 = 10
num_plots_0 = a.shape[0]

line_style = '-'

plt.figure()

### plot along a line ###
# plt.subplot(2,1,2)
# plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,num_plots_1)])
# for ii in range(num_plots_1):
#     plt.plot(a[:,ii],line_style,label=str(ii))
#     plt.title('along line x=y at time :') 
# plt.legend()


plt.subplot(2,1,1) 
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,num_plots_0)])
for ii in range(num_plots_0):
    plt.plot(a[ii,:],line_style,label=str(ii)) 
    plt.title('time history of cell : ') 
# plt.legend()


### plot the spectrum ###

plt.subplot(2,1,2)
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,num_plots_0)])
for ii in range(num_plots_0):

    print ii
    print a.shape
    data = np.reshape(a[ii,:], (a.shape[1],1))
    freq, ffty = ComputeSpectrum( 80000., data, N=128 ) 
    plt.loglog(freq,ffty,line_style,label=str(ii))
    plt.title('frequency response of cell : :') 
# plt.legend()


plt.show()
