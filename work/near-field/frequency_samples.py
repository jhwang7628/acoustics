#!/usr/bin/env python

import sys
import os
from subprocess import call

def ConfigString(f, outputDirectory): 

    s = """<?xml version="1.0" ?>
<impulse_response>
    <mesh file="/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj" scale="1.0" />

    <solver
        stop_time="0.01"
        distancefield="/home/jui-hsien/code/acoustics/work/meshes/small_ball/small_ball.obj.1m.dist"
        fieldresolution="100"
        gridresolution="100"
        cellsize="0.01" 
        timestepfrequency="176400"
        substeps="4"

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

</impulse_response> """  %(f,outputDirectory,'%s')

    return s

def main(): 

    frequencySamples = range(500,22100,500)
    print frequencySamples
    
    
    for f in frequencySamples: 
    
        outputDirectory = 'ball_experiment/%uHz' %(f)
    
        if os.path.isdir(outputDirectory): 
            print 'directory %s exists, skipping' %(outputDirectory)
            continue
        
        print '\n\n=== sampling %s ===' %(outputDirectory) 

        print ConfigString(f, outputDirectory)


        of = open('temp_config.xml','w')
        of.write(ConfigString(f, outputDirectory)) 
        of.close()
        break


if __name__ == '__main__': 
    main()


