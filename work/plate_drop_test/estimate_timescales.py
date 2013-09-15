#!/usr/bin/python
#
# Estimates impact timescales using Hertz contact theory

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

impulseFile = 'impulses.txt';
forceScalePrefix = 'forceScales';

cmd = '%s/estimate-timescales %s %s' \
       % ( binPath, impulseFile, forceScalePrefix );

print( "=== Estimating Hertz contact time scales" );
os.system(cmd);
print( "=== Done" );
