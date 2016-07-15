#!/usr/bin/python
# create tet mesh from obj, the parameters work well for 'plate' model

import sys
import os

if len(sys.argv) != 3:
    print '**Usage: %s <input_obj_file> <output_tet_file>' %(sys.argv[0])
    sys.exit()

# parse input arguments
objFile = sys.argv[1]
tetFile = sys.argv[2]

installPath = '/home/jui-hsien/code/acoustics';
binPath = '%s/build_release/bin' %(installPath)
redirect = '/dev/stdout'

# define parameters and call isostuffer
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
