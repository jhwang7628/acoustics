#!/usr/bin/env python
import numpy as np 
import matplotlib.pyplot as plt
import struct,glob,time
##
datadir = 'data'
filenames = sorted(glob.glob('%s/*.dat' %(datadir)))
##
np.set_printoptions(precision=3)
fig = plt.figure()
for filename in filenames: 
    print '------------ %s --------------' %(filename)
    ifs = open(filename, 'rb')
    size = struct.unpack('ii', ifs.read(8))
    data = np.zeros((size[0],size[1]))
    for col in range(size[1]): 
        for row in range(size[0]): 
            data[row, col] = struct.unpack('d', ifs.read(8))[0]
    if filename == filenames[0]:
        p = plt.imshow(data, interpolation='none', clim=(-10.,10.))
    else: 
        p.set_data(data)
    plt.pause(0.01)


