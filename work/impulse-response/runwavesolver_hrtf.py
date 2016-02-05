#!/usr/bin/env python 

import subprocess
import os
import numpy as np
from runwavesolver import RunWavesolver

def CreateHRTFConfigString(left_right, source_x, source_y, source_z): 

    configString = """<?xml version="1.0" ?>
    <impulse_response>
        <mesh file="/home/jui-hsien/code/acoustics/work/impulse-response/geometry/fluent5.msh.obj" scale="1.0" />
    
        <solver
            stop_time="0.005"
            distancefield="/home/jui-hsien/code/acoustics/work/impulse-response/geometry/fluent5.msh.dist"
            fieldresolution="250"
            gridresolution="250"
            gridscale="1.0800000"
            timestepfrequency="640000"
            substeps="8"
    
            output_pattern="data_Gaussian_space_depend_time_1over160000_res250_simulation/head_%s_%s"
    
            source_position_x="%f" 
            source_position_y="%f" 
            source_position_z="%f" 
            />
    
    </impulse_response>
    """ %(left_right, '%s', source_x, source_y, source_z) 

    return configString

def main(): 

    fileID = 0
    while os.path.isfile('head_left_%.5u.xml' %(fileID)) or os.path.isfile('head_right_%.5u.xml' %(fileID)): 
        fileID+=1
   
    for ii in range(2): 
        source_x = 0.0
        source_y = 0.0
        source_z = 0.0
        if ii==0: 
            left_right = 'left' 
            source_z = 0.08
        else: 
            left_right = 'right'
            source_z = -0.08
        configString = CreateHRTFConfigString(left_right, source_x, source_y, source_z) 
        configFile = 'head_%s_%.5u.xml' %(left_right, fileID)
        outFile = open(configFile,'w')
        outFile.write(configString)
        outFile.close()

        RunWavesolver(configFile)
        # subprocess.call('rm %s' %(configFile), shell=True)

if __name__=="__main__": 
    main()
