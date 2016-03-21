#!/usr/bin/env python 

import my_io as io 
import sys

if (len(sys.argv)!=3): 
    print '**Usage: %s <vertex_position_file> <query_cell_index>' %(sys.argv[0])
    sys.exit()

a = io.readMatrixXdBinary(sys.argv[1])
print 'query cell %u has position' %(int(sys.argv[2])), a[int(sys.argv[2]),:]
