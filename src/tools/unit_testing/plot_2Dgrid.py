#!/usr/bin/env python
import numpy as np 
import matplotlib.pyplot as plt
import struct,glob,time,sys
if len(sys.argv) != 2: 
    print '**Usage: %s <datadir>' %(sys.argv[0])
    sys.exit()
##
datadir = sys.argv[1]
filenames = sorted(glob.glob('%s/*.dat' %(datadir)))
stop = None
##
np.set_printoptions(precision=3)
fig = plt.figure()
listened = []
count = 0
for filename in filenames: 
    print '------------ %s --------------' %(filename)
    ifs = open(filename, 'rb')
    size = struct.unpack('ii', ifs.read(8))
    data = np.zeros((size[0],size[1]))
    for col in range(size[1]): 
        for row in range(size[0]): 
            data[row, col] = struct.unpack('d', ifs.read(8))[0]
    if filename == filenames[0]:
        p = plt.imshow(data, interpolation='none', clim=(-50,50))
    else: 
        p.set_data(data)
    fetch_pt = (int(size[0]/2+70), int(size[1]/2))
    print fetch_pt
    ldata = data[fetch_pt[0], fetch_pt[1]]
    # ldata = data[int(size[0]/2), int(size[1]/2)]
    listened.append(ldata)
    plt.pause(0.01)
    # plt.savefig('frames/test2-%.4d.jpg' %(count))
    count += 1
    if stop and count == stop:
        break
plt.figure()
plt.plot(listened)
plt.show()

# save
ofs = open('ldata-%s.txt' %(datadir), 'w') 
for l in listened: 
    ofs.write('%.16f\n' %(l))
ofs.close()
