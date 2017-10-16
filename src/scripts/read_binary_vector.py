#!/usr/bin/env python
import binary_io as io
import sys

if len(sys.argv) != 2: 
    print '**Usage: %s <vector_filename>' %(sys.argv[0]) 
    sys.exit()

filename = sys.argv[1]
data = io.readVectorBinary(filename)
print 'Vector read has dimension: ', data.shape
print data

