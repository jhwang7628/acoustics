#!/usr/bin/env python 
import sys
import glob
from subprocess import call

## NOTE: only execute this in the data folder
if len(sys.argv) != 3: 
    print '**Usage: %s <old_prefix> <new_prefix>'
    sys.exit(1)

oldPrefix = sys.argv[1]
newPrefix = sys.argv[2]
filenames = glob.glob('%s*' %(oldPrefix))
for f in filenames: 
    cmd = 'mv %s %s' %(f, f.replace(oldPrefix, newPrefix))
    print cmd
    call(cmd, shell=True)
