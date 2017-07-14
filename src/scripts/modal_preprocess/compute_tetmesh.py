#!/usr/bin/env python
import sys
import os

################################################################################
# This code create tet mesh from obj, the parameters work well for 'plate' model
################################################################################
if len(sys.argv) < 3:
    print '**Usage: %s <input_obj_file> <output_tet_file> [isoResolution=6]' %(sys.argv[0])
    sys.exit()

## User defined settings
installPath = '/home/jui-hsien/code/acoustics';
binPath = '%s/build_release/bin' %(installPath)
redirect = '/dev/stdout'

isoNlevel = 5; # number of levels
isoMargin = 1; # margin of the bounding box
isoAlpha = 0.25; # alpha long value in isostuff algorithm
isoBeta = 0.42978; # alpha short value in isostuff algorithm

## Automatic
objFile = sys.argv[1]
tetFile = sys.argv[2]
if len(sys.argv) == 4:
    isoResolution = int(sys.argv[3]) # bounding box resolution
else:
    isoResolution = 6;

## Pipeline starts
print '################################################################################'
print '## Generate Tet mesh'
print '################################################################################'
cmd = '%s/isostuffer -R %d -L %d -M %d -a %f -b %f %s %s' \
       % ( binPath, isoResolution, isoNlevel, isoMargin, isoAlpha, isoBeta, \
           objFile, tetFile );
print cmd;
print '........................................'
os.system( '%s > %s' % ( cmd, redirect ) );
print '\n\n'
