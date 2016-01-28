#!/usr/bin/env python 
import numpy as np
import matplotlib.pyplot as plt



a = np.loadtxt('tmp.txt') 



plt.figure()
plt.plot(a[:,0],a[:,1])
plt.show(block=True)

