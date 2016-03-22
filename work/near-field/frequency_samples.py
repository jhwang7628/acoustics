#!/usr/bin/env python

import sys
import os
from subprocess import call
import runwavesolver
import time
import numpy as np

def ConfigString(f, outputDirectory, endtime=0.01, timestepfrequency=176400, substeps=4): 

    s = """<?xml version="1.0" ?>
<impulse_response>
    <mesh file="/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj" scale="1.0" />

    <solver
        stop_time="%f"
        distancefield="/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj.1m.dist"
        fieldresolution="150"
        gridresolution="150"
        cellsize="0.0066666666666666667" 
        timestepfrequency="%u"
        substeps="%u"

        f="%u"
        output_pattern="%s/test_%s"

        use_mesh="1"
        cornell_box_boundary_condition="0"

        /> 

    <volumetric_source> 
        <source
        source_width_time="0.00005"
        source_position_x="0.0" 
        source_position_y="0.0" 
        source_position_z="0.0" 
        />
    </volumetric_source> 

</impulse_response> """  %(endtime, int(timestepfrequency), int(substeps), f,outputDirectory,'%s')

    return s

def main(): 

    # frequencySamples = range(500,22100,500)

    frequencySamples = np.logspace(np.log10(100),np.log10(22100),30)
    frequencySamples = frequencySamples.astype(int)
    print 'frequency samples: ', frequencySamples
    
    
    for f in frequencySamples: 
    
        outputDirectory = 'data_ball_experiment/res150/%uHz' %(f)
   

        # debug f range
        if os.path.isdir(outputDirectory): # and f < 10000: 
            print 'directory %s exists, skipping' %(outputDirectory)
            continue

        
        print 'sampling %s' %(outputDirectory) 


        configFile = 'solver_config_%u.xml' %(f)
        of = open(configFile, 'w')

        # debug 
        timestepfrequency=max(f*40*4, 176400)
        endtime=min(0.0007+1./f*10.,0.01)
        of.write(ConfigString(f, outputDirectory, endtime=endtime, timestepfrequency=timestepfrequency, substeps=4 ))  # make sure all sine has at least 20 points
        of.close()


        # debug 
        # call('rm %s/*' %(outputDirectory), shell=True)

        t0 = time.time()
        redirect = '%s/log.runwavesolver' %(outputDirectory)
        runwavesolver_strings = ['runwavesolver', configFile, '1']
        runwavesolver.main(runwavesolver_strings, redirect)
        t1 = time.time()

        print ' %f seconds for wavesolver run' %(t1-t0)

        print ' removing xml config'
        call('rm %s' %(configFile), shell=True)


if __name__ == '__main__': 
    main()


