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

# Put the frames in ./frames
frameDir = 'frames';
os.system('mkdir -p %s' % ( frameDir ) );

configFile = 'default.cfg';
displaceFile = 'displace.bin';

# The rigid simulation ran at 1000Hz.  We show every 10th time step, which
# will result in a frame rate of 100 FPS.
frameFrequency = 10;

# Bring up the replayer
cmd = '%s/replayer -f %s -d %s -r %d -o %s/frame' \
       % ( binPath, configFile, displaceFile, frameFrequency, frameDir );

print( "=== Runing the rigid-body replayer" );
print( "=== To generate frames:" );
print( "===     1) Adjust the camera view as desired" );
print( "===     2) Press 'P' to enable frame snapshots" );
print( "===     3) Press space to start writing frames" );
print( "" );

os.system(cmd);
