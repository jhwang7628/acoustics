#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt
import sys
filenames = sys.argv[1:]
plt.figure(figsize=(12,6))
data = dict()
for f in filenames: 
    key = f.split('_')[-1].split('.')[0]
    data[key] = np.loadtxt(f) 
    plt.plot(data[key], label='%s' %(key), linewidth=1.5)
plt.legend(loc=2)

ref = 'baseline'
ref_val = max(data['baseline'])
save_name = 'plot'
title = 'Error percentage: '
for d in data: 
    save_name += '-%s' %(d)
    if d != ref: 
        err = (data[d][-1]-data[ref][-1])/ref_val
        print '\'%s\' end error compared to baseline \'%s\': %.8f' %(d, ref, err)
        title += '%s=%.3f; ' %(d, err*100.)
save_name += '.pdf'
plt.xlabel('time step')
plt.ylabel('pressure at cell near padded boundary')
plt.title(title)
plt.savefig(save_name)
plt.show()
