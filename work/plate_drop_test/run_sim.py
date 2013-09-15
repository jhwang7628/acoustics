#!/usr/bin/python
#
# Script which runs the initial rigid-body simulation

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

cmd = '%s/rigidsim default.cfg' % ( binPath );
print( "=== Running the rigid-body simulator" );
os.system(cmd);
print( "=== Done" );
