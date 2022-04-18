#!/usr/bin/env python
import sys,os
import matplotlib.pyplot as plt
from wavesolver_results import *
import scipy

if len(sys.argv) != 2: 
    print('**Usage: %s <data_folder>' %sys.argv[0])
    sys.exit()

##
folder = sys.argv[1]
##
results = Wavesolver_Results()
results.Set_Folder(folder)
all_data = results.Read_All_Audio()
N_points = all_data.shape[1]
N_steps  = all_data.shape[0] 

plt.figure()
for ii in range(N_points): 
    plt.plot(all_data[:,ii])
plt.show()
