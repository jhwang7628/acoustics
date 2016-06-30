#!/usr/bin/env python 

import my_io as io 
import sys

if (len(sys.argv) != 2): 
    print '**Usage: %s <binary_dat_file>' %(sys.argv[0])
    sys.exit()

a = io.readMatrixXdBinary(sys.argv[1])
print a
