#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt 

a = np.linspace( 0.07, 0.21899, 100 )
aa = np.column_stack((np.zeros(a.size),np.zeros(a.size),a))
np.savetxt('head_listening.txt', aa)
