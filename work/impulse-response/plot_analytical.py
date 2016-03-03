#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt 
import sys 


if len(sys.argv) < 2: 
    print '**Usage: %s <analytical_data>' %(sys.argv[0])
    sys.exit()


plt.figure()
line_styles = ['-', 'o']
count=0
for filename in sys.argv[1:]:

    l = np.loadtxt('head_listening.txt') 
    # l = np.loadtxt(filename) 
    
    analyticalData = np.loadtxt(filename)[:,:]
    
    colormap = plt.cm.jet
    line_style = line_styles[count]
    
    N = analyticalData[:,:50].shape[1]
    print N
    
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0,0.9,N)])
    for ii in range(N):
        plt.plot(l[:,-1], analyticalData[:,ii],line_style,label=str(ii))
    plt.xlabel('x')
    plt.ylabel('pressure')
    count+=1

plt.show()


