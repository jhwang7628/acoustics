#!/usr/bin/python
#
# Runs (some of) the steps needed to precompute modal sound data for an object

import sys
import os

objName = 'plate';
objFile = 'plate.obj'
tetFile = 'plate.tet';
dataPrefix = 'plateRigid';

# By default, don't product output, just write the commands
redirect = '/dev/null';
#redirect = '/dev/stdout';

youngsModulus = 72000000000.0;
poissonRatio = 0.19;
density = 2300.0;

# More than enough for the plate
numEigs = 100;


installPath = '/home/jchadwick/research/acoustics';
if ( len(sys.argv) >= 2 ):
  installPath = sys.argv[1];

binPath = '%s/gcc-build/src' % ( installPath );

matlabPath = '%s/src/matlab' % ( installPath );

# Step one - generate the tet mesh
isoResolution = 7;
isoNlevel = 3;
isoMargin = 7;
isoAlpha = 0.25;
isoBeta = 0.42978;

cmd = '%s/isostuffer -R %d -L %d -M %d -a %f -b %f %s %s' \
       % ( binPath, isoResolution, isoNlevel, isoMargin, isoAlpha, isoBeta, \
           objFile, tetFile );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );

# Step two - generate the mass and stiffness matrix
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

# Step four - generate input files for the Helmholtz radiation problems
os.system( 'mkdir -p fbem_in' );
cmd = '%s/fbem_input_gen -m %s.modes -d %f -t %s -o fbem_in/input-%%d.txt' \
       % ( binPath, objName, density, tetFile );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );

# Step five - write the object center of mass, rigid inertia matrix, and
# material parameters to files that will be needed in other parts of the
# pipeline
cmd = '%s/write-mass-center %s %s' % ( binPath, tetFile, dataPrefix );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );

cmd = '%s/write-inertia %s %s' % ( binPath, tetFile, dataPrefix );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );

# This part assumes that octave is installed.  This command could be
# changed somewhat to work with Matlab
cmd = ("octave --eval \"addpath('%s'); "
       "make_material_parameters(%f, %f, '%s_materialProperties.vector');\"") \
       % ( matlabPath, youngsModulus, poissonRatio, dataPrefix );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );

# Step six - wrap everything up in a tarball so that it can be conveniently
# moved to the hydra cluster.
#
# As a final step, move this file to /hydra/S1/impact_snd/fastbem/models and
# untar it there (e.g., tar xzvf plate.tar.gz).  This should create a 'plate'
# directory.  Inside this directory, you can run the following command:
#     ./../../scripts/fbem_solve.pl . <modeID>
# This will solve the boundary element problem for the given mode ID, placing
# the results in fbem_ret
# 
cmd = ("tar --transform 's,^,%s/,' -czvf %s.tar.gz fbem_in/ "
       "freqs.txt %s.geo.txt") % ( objName, objName, objName );
print cmd;
os.system( '%s > %s' % ( cmd, redirect ) );
