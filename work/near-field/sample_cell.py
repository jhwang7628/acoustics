#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np 
import my_io as io
import glob 
import sys


if len(sys.argv) < 2: 
    print '**Usage: %s <folders..>' %(sys.argv[0])
    sys.exit()


N_dir=len(sys.argv[1:])

frequency = np.zeros(len(sys.argv[1:]))
magnitude = np.zeros(len(sys.argv[1:]))
analytical2 = np.ones(len(frequency))

plt.figure()

count2 = 0
for d in sorted(sys.argv[1:]): 

    print 'processing %s' %(d) 

    tokens = d.split('_')
    print tokens
    for t in tokens: 
        if t.find('Hz') >= 0: 
            print t
            f = float(t[:-2])

    if f<10000: 
        continue 

    frequency[count2] = f

    filenames = glob.glob('%s/*pressure*.dat' %(d))
    filenames = sorted(filenames)
 
    ## debug for the ball 100^3
    ##        509550   509450
    ##  -----------------
    ##  |  solid| fluid |
    ##  |   o   |   o   |
    ##  |504550 | 504450|
    ##  -----------------
    ## velocity interfacial cell: 509550 -> 0.005 -0.05 0.005
    ##                            504450 -> pressure cell right next to it
    cell_id = 505850 # -0.005,-0.285,-0.005
    # cell_id = 255050 # for 100
    # cell_id = 866325 # for 150
    # cell_id = 502150 # -0.005,-0.285,-0.005
    # cell_id = 4008500 # 0.0025, -0.2875, 0.0025
    signal=np.zeros(len(filenames))
    analytical=np.zeros(len(filenames))


    vertex_position_filename = '%s/test_vertex_position.dat' %(d) 
    pos = io.readMatrixXdBinary(vertex_position_filename, verbose=0)[cell_id,:] 
    print 'sampled cell has position: ', pos.transpose()


    
    count = 0
    for f in filenames: 
        tmp = io.readMatrixXdBinary(f,verbose=0)
        signal[count] = tmp[cell_id]

        # for plane 
        # analytical[count] = 343.

        # for sphere
        analytical2[count2] = 0.05**2/np.linalg.norm(pos)/np.sqrt(1+(2.*np.pi*frequency[count2]*0.05/343.)**2)
        # analytical2[count2] = 0.05**2/np.sqrt(0.245**2+2.*0.005**2)/np.sqrt(1+(2.*np.pi*frequency[count2]*0.05/343.)**2)
        # analytical[count] = 0.05**2/np.sqrt(0.365**2*3)
        # analytical[count] = 0.05**2/np.sqrt(0.2875**2+2.*0.0025**2)
        count += 1

    magnitude[count2] = max(signal)

    t = range(len(filenames))
    plt.plot(t, signal, label=d)
    count2 += 1

plt.legend()
# plt.plot(t, -analytical, 'k--', label='analytical')

plt.figure()
plt.subplot(2,1,1)
plt.plot(frequency, magnitude, 'o', label='numerical') 
plt.xlabel('frequency (Hz)') 
plt.ylabel('pressure magnitude (Pa)') 
# print analytical
plt.plot(frequency, analytical2, 'xk', label='analytical')
plt.xlabel('frequency (Hz)') 
plt.ylabel('pressure magnitude (Pa)') 
plt.grid()
plt.legend()


# plt.plot(frequency, analytical, 'k--', label='analytical')
plt.subplot(2,1,2)
plt.loglog(frequency, magnitude, 'o', label='numerical log') 
plt.loglog(frequency, analytical2, 'xk', label='analytical log')
plt.xlabel('frequency (Hz)') 
plt.ylabel('pressure magnitude (Pa)') 
plt.grid()
plt.legend()

plt.figure()
plt.plot(frequency, (magnitude-analytical2)/analytical2, 'o', label='relative error')
plt.xlabel('frequency (Hz)') 
plt.ylabel('numerical - analytical / analytical') 
plt.grid()

plt.legend()

plt.show()



