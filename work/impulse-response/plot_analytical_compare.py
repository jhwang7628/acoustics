#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt 
import sys 


if len(sys.argv) != 3: 
    print '**Usage: %s <computed_data> <analytical_data>' %(sys.argv[0])
    sys.exit()

l = np.loadtxt('../head_listening.txt') 

chop=15
chopend=1

computedData = np.loadtxt(sys.argv[1])[chop:-chopend,:]
analyticalData = np.loadtxt(sys.argv[2])[chop:-chopend,:]

dt_computed   = 1./160000.
endTime_computed    = dt_computed*float(computedData.shape[1])
endTime_analytical  = endTime_computed 
dt_analytical = endTime_computed/float(analyticalData.shape[1])

num_plots_1_computed = computedData.shape[0]
num_plots_0_computed = computedData.shape[1]

num_plots_1_analytical = analyticalData.shape[0]
num_plots_0_analytical = analyticalData.shape[1]

colormap = plt.cm.jet
line_style = '-'

plt.figure()

t_computed   = np.linspace(0.0,endTime_computed,computedData.shape[1])
t_analytical = np.linspace(0.0,endTime_analytical,analyticalData.shape[1])

x = l[chop:-chopend]

plt.subplot(2,2,1)
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,num_plots_0_computed)])
for ii in range(num_plots_0_computed):
    plt.plot(x[:,2], computedData[:,ii],line_style,label=str(ii))
plt.xlabel('x')
plt.ylabel('pressure')
plt.title('Wave solver: along line x=y at given time') 
# print plt.axis()
# plt.legend()

plt.subplot(2,2,2) 
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,num_plots_1_computed)])
for ii in range(num_plots_1_computed):
    plt.plot(t_computed, computedData[ii,:],line_style,label=str(ii)) 
plt.xlabel('t')
plt.ylabel('pressure')
plt.title('Wave solver: time history of cell ') 
# print plt.axis()
# plt.legend()

plt.subplot(2,2,3)
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,num_plots_0_analytical)])
for ii in range(num_plots_0_analytical):
    plt.plot(x[:,2], analyticalData[:,ii],line_style,label=str(ii))
plt.xlabel('x')
plt.ylabel('pressure')
plt.title('Analytical: along line x=y at given time') 
# plt.legend()

plt.subplot(2,2,4) 
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,num_plots_1_analytical)])
for ii in range(num_plots_1_analytical):
    plt.plot(t_analytical, analyticalData[ii,:],line_style,label=str(ii)) 
plt.xlabel('t')
plt.ylabel('pressure')
plt.title('Analytical: time history of cell ') 
# plt.legend()
 
plt.show()





