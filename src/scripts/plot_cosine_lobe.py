#!/usr/bin/env python
import matplotlib.pyplot as plt 
import numpy as np 
import sys
from my_io import readMatrixXdBinary

if len(sys.argv) != 4: 
    print '**Usage: %s <filename_GC> <filename_RAS> <listening_position_file>' %(sys.argv[0])
    sys.exit()

filename_GC = sys.argv[1]
filename_RAS= sys.argv[2]
posFile = sys.argv[3]

data_GC = readMatrixXdBinary(filename_GC)
data_RAS = readMatrixXdBinary(filename_RAS)
positions = readMatrixXdBinary(posFile)
N_points = data_GC.shape[0]

theta = np.arctan2(positions[:,1], positions[:,0])
lobePositions_GC = np.zeros((N_points, 2))
lobePositions_RAS = np.zeros((N_points, 2))
theoreticalLobePositions = np.zeros((N_points, 2))

lobePositions_GC[:, 0] = np.multiply(np.square(data_GC[:,0]), np.cos(theta))
lobePositions_GC[:, 1] = np.multiply(np.square(data_GC[:,0]), np.sin(theta))
lobePositions_RAS[:, 0] = np.multiply(np.square(data_RAS[:,0]), np.cos(theta))
lobePositions_RAS[:, 1] = np.multiply(np.square(data_RAS[:,0]), np.sin(theta))
theoreticalLobePositions[:, 0] = np.multiply(np.square(np.cos(theta)), np.cos(theta))
theoreticalLobePositions[:, 1] = np.multiply(np.square(np.cos(theta)), np.sin(theta))
maxLobe_GC = lobePositions_GC.max()
maxLobe_RAS = lobePositions_RAS.max()
print 'maxLobe_GC = %f; maxLobe_RAS = %f.' %(maxLobe_GC, maxLobe_RAS)
# lobePositions_GC /= maxLobe_GC
# lobePositions_RAS /= maxLobe_RAS

fig = plt.figure(figsize=[10, 10])
lw = 2.
plt.plot(theoreticalLobePositions[:,0], theoreticalLobePositions[:,1], 'k--', linewidth=lw)
plt.plot(lobePositions_GC[:,0], lobePositions_GC[:,1], 'b', linewidth=lw, label=filename_GC)
plt.plot(lobePositions_RAS[:,0], lobePositions_RAS[:,1], 'r', linewidth=lw, label=filename_RAS)
plt.title('Directivity of AN radiation (dipole). Resolution: 0.5cm')
plt.xlabel('x')
plt.xlabel('y')
plt.legend(loc=0)
plt.axis('equal')
# plt.axis([-1, 1, -1, 1])
plt.grid() 
plt.show()
