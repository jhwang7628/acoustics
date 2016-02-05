import argparse
import os
import sys
from subprocess import call
import xml.etree.ElementTree as ET


def RunWavesolver(configFile):

    # configFile = args[0]

    print configFile
    tree = ET.parse(configFile)
    root = tree.getroot()

    pattern = root.find('solver').get('output_pattern')
    outdir,outname = os.path.split(pattern)

    if (not os.path.isdir(outdir)): 
        warningMessage = '**WARNING** requested directory %s does not exist. create one and continue? [y/N] ' %(outdir)
        createDir = raw_input(warningMessage) 
        if (createDir == 'y'): 
            call('mkdir -p %s' %(outdir), shell=True)
        else: 
            sys.exit()

    call('cp %s %s' %(configFile, outdir), shell=True)

    os.environ["PATH"] += os.pathsep + ("/home/jui-hsien/code/acoustics/build/bin")
    os.system("precompute-impulse-response-text %s"%(configFile))

# simple wrapper that deals with the directory check, etc
if __name__ == '__main__':
    if (len(sys.argv)!=2): 
        print '**Usage %s <impulse_response_config_file>' %(sys.argv[0])
        sys.exit()

    RunWavesolver(sys.argv[1])
