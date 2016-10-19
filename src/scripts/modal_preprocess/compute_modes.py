#!/usr/bin/env python
import sys
import os
################################################################################
## This code create tet mesh from obj
################################################################################
if len(sys.argv) != 2: 
    print '**Usage: %s <obj_prefix>' %(sys.argv[0])
    sys.exit()

## User defined settings
installPath = '/home/jui-hsien/code/acoustics'
binPath = '%s/build_release/bin' %(installPath)
redirect = '/dev/stdout';
                     
youngsModulus = 72000000000.0;
poissonRatio = 0.19;
density = 2300.0;
           
numEigs = 100;

## Automatic 
objName = sys.argv[1];
tetFile = '%s.tet' %(objName)

## Pipeline starts
# Generate the mass and stiffness matrix
cmd = '%s/elasticity_solver %s %f %f' %(binPath, tetFile, youngsModulus, poissonRatio);
print cmd;
os.system('%s > %s' %(cmd, redirect));

# This step also generates <objName>.tet.geo.txt, but we want to rename this
os.system('mv -f %s.geo.txt %s.geo.txt' %(tetFile, objName));

# Linear modal analysis
# Using a threshold of 1.0 here seems to work
cmd = ("%s/arpack-eigensolver -n %d -t 1.0 -s %s_stiffness.mat -m %s_mass.mat -o %s.modes -v") %(binPath, numEigs, tetFile, tetFile, objName);
print cmd;
os.system('%s > %s' %(cmd, redirect));
