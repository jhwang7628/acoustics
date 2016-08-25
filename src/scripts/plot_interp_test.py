#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt

rawFile = 'rawSignal.txt'
interpFile = 'interpolatedSignal.txt'

print 'read files'
raw = np.loadtxt(rawFile)
interp = np.loadtxt(interpFile) 

print 'plotting'
plt.figure()
plt.plot(raw[:,0], raw[:,1], 'bo')
plt.plot(interp[:,0], interp[:,1], 'rp')
plt.show()

