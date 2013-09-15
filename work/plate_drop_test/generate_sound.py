#!/usr/bin/python
#
# Generates sound from rigid-body simulation data

import sys
import os

installPath = '/home/jchadwick/research/acoustics';
if ( len(sys.argv) >= 2 ):
  installPath = sys.argv[1];

# Optionally specify where things are built
if ( len(sys.argv) >= 3 ):
  binPath = '%s/%s/src' % ( installPath, sys.argv[2] );
else :
  binPath = '%s/gcc-build/src' % ( installPath );

matlabPath = '%s/src/matlab' % ( installPath );

configFile = 'example.cfg';

cmd = '%s/sound-generator %s' % ( binPath, configFile );

print( "=== Running modal sound generator" );
os.system(cmd);
print( "=== Done" );

# Use octave to turn the output in to a wave file
cmd = "octave --eval \"addpath('%s'); make_sound();\"" % ( matlabPath );
os.system(cmd);
