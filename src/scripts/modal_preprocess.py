#!/usr/bin/python
# Runs (some of) the steps needed to precompute modal sound data for an object
import sys
import os

from get_paths import *

################################################################################
## set this block for the current run
################################################################################
## required
# fs and install paths
objName = 'big_plane';
installPath = projectPath()
binPath = '%s/build_release/bin' %(installPath)

## run parameters
# material properties
youngsModulus = 72000000000.0;
poissonRatio = 0.19;
density = 2300.0;

# tet mesh properties
# isoResolution = 7;
isoResolution = 5;
isoNlevel = 3;
isoMargin = 7;
isoAlpha = 0.25;
isoBeta = 0.42978;

# modal properties
numEigs = 100;

## optional
matlabPath = '%s/src/matlab' % ( installPath );
redirect = '/dev/stdout';
################################################################################


################################################################################
## Below are three stages needed for modal analysis. 
## Commented out sections if not needed.
################################################################################
objFile = '%s.obj' %(objName)
tetFile = '%s.tet' %(objName)
dataPrefix = '%sRigid' %(objName)

# Step one - generate the tet mesh
cmd = '%s/isostuffer -R %d -L %d -M %d -a %f -b %f %s %s' \
       % ( binPath, isoResolution, isoNlevel, isoMargin, isoAlpha, isoBeta, \
           objFile, tetFile );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );

# # Step two - generate the mass and stiffness matrix
cmd = '%s/elasticity_solver %s %f %f' \
       % ( binPath, tetFile, youngsModulus, poissonRatio );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );

# This step also generates plate.tet.geo.txt, but we want to rename this
os.system( 'mv -f %s.geo.txt %s.geo.txt' % ( tetFile, objName ) );

# Step three - linear modal analysis
#
# Using a threshold of 1.0 here seems to work
cmd = ("%s/arpack-eigensolver -n %d -t 1.0 -s %s_stiffness.mat -m %s_mass.mat "
       "-o %s.modes -v") % ( binPath, numEigs, tetFile, tetFile, objName );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );
