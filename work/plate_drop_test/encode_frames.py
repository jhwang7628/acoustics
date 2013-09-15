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

# The frames are in ./frames
frameDir = 'frames';

encodeScript = '%s/src/scripts/encode_png_sequence.py' % ( installPath );

# We ran a rigid-body simulation at 1000Hz, then dumped out every 10th time
# step, resulting in a rate of 100 FPS
frameRate = 100;

videoFile = 'test.mpg';

imageXres = 1024;
imageYres = 576;

# Make a movie!
cmd = 'python %s %s/frame %s %d %d %d' \
       % ( encodeScript, frameDir, videoFile, imageXres, imageYres, frameRate );

os.system(cmd);
