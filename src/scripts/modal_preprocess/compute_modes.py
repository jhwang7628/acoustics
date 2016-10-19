#!/usr/bin/env python
import sys
import os
################################################################################
## This code create tet mesh from obj
################################################################################
if len(sys.argv) != 3: 
    print '**Usage: %s <obj_prefix> <num_eigenvalues>' %(sys.argv[0])
    sys.exit()

## User defined settings
installPath = '/home/jui-hsien/code/acoustics'
binPath = '%s/build_release/bin' %(installPath)
redirect = '/dev/stdout';
                     
youngsModulus = 72000000000.0;
poissonRatio = 0.19;
density = 2300.0;
           
## Automatic 
objName = sys.argv[1];
numEigs = int(sys.argv[2]);
tetFile = '%s.tet' %(objName)
bin_elasticity='%s/elasticity_solver' %(binPath)
bin_arpack_eigensolver='%s/arpack-eigensolver' %(binPath)
if not os.path.isfile(bin_elasticity) or not os.path.isfile(bin_arpack_eigensolver): 
    print '**ERROR** some binaries required are not found in path %s' %(binPath)

## Pipeline starts
print '################################################################################'
print '## Generate mass and stiff matrix'
print '################################################################################'
cmd = '%s %s %f %f' %(bin_elasticity, tetFile, youngsModulus, poissonRatio);
print cmd;
print '........................................'
os.system('%s > %s' %(cmd, redirect));
print '\n\n'

# This step also generates <objName>.tet.geo.txt, but we want to rename this
os.system('mv -f %s.geo.txt %s.geo.txt' %(tetFile, objName));

print '################################################################################'
print '## Linear modal analysis'
print '################################################################################'
# Using a threshold of 1.0 here seems to work
cmd = ("%s -n %d -t 1.0 -s %s_stiffness.mat -m %s_mass.mat -o %s.modes -v") %(bin_arpack_eigensolver, numEigs, tetFile, tetFile, objName);
print cmd;
print '........................................'
os.system('%s > %s' %(cmd, redirect));
print '\n\n'
