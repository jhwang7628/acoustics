#!/usr/bin/env python

################################################################################
## This script reads the file stored by MAC_Grid::ToFile_GhostCellCoupledMatrix()
################################################################################

import sys
import numpy as np 
import matplotlib.pyplot as plt

if len(sys.argv) != 2: 
    print '**Usage: %s <file_name>' %(sys.argv[0])
    sys.exit()

inFile = open(sys.argv[1])
N_ghostCells = int(inFile.readline())
inFile.readline()
mat = np.zeros((N_ghostCells, N_ghostCells)) 
plt.figure()

print 'N = %u' %(N_ghostCells)
print 'Reading and parsing'
for row in range(N_ghostCells): 
    if row % 500 == 0:
        print 'row = ', row
    Nnz = int(inFile.readline()) 
    for nz in range(Nnz): 
        line = inFile.readline().split()
        col = int(line[0]) 
        val = float(line[1]) 
        plt.plot(row, col, 'bo')
        mat[row, col] = val

plt.axis('equal')
plt.axis([0, N_ghostCells, 0, N_ghostCells])
plt.show()

# print 'Drawing'
# plt.figure()
# plt.pcolor(mat[::20, ::20], cmap='RdBu', vmin=-1, vmax=1)
# plt.colorbar()
# plt.show()
