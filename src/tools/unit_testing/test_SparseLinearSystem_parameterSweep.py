#!/usr/bin/env python 

import numpy as np 
from subprocess import call


testBinary = '/home/jui-hsien/code/acoustics/build_test/bin/unit_testing_sparseLinearSystemSolver'

numberTest = 500
testDense = 0


## first test on fixed dim with varying nnz
dim = 1000
nnz = np.logspace(np.log10(10), np.log10(999999),num=numberTest)

for n in nnz: 
    n = int(n)
    filename = 'timing_fixeddim_biCGStab.txt' 
    call('%s %u %u %u %s' %(testBinary, dim, n, testDense, filename),shell=True) 



## second test on fixed nnz with varying dim
# dim = np.logspace(np.log10(25),np.log10(1000000),num=numberTest)
# nnz = 500
# 
# for d in dim: 
#     d = int(d)
#     filename = 'timing_fixednnz_biCGStab.txt' 
#     call('%s %u %u %u %s' %(testBinary, d, nnz, testDense, filename),shell=True) 





