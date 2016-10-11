#!/usr/bin/env python
import numpy as np
import glob
from subprocess import call

def Render(cmd, templateFile, tmpFilename, outName, objectXMLString): 
    # first write a tmp file
    of = open(tmpFilename, 'w')
    of.write(objectXMLString)
    of.close()
    cmd = cmd %(tmpFilename, outName)
    print cmd
    call(cmd, shell=True)

def CreateObjectXML(x, y, z): 
    s  = """<scene version="0.5.0">"""
    s += """
    <shape type="obj">
        <string name="filename" value="assets/sphere_d10cm.obj"/>
        <transform name="toWorld">
            <translate x="%f" y="%f" z="%f"/>
        </transform>
        <ref id="sphere"/>
    </shape>\n""" %(x, y, z)
    s += """</scene>"""
    return s

def ExtendLastFrame(filenames_s, append=20): 
    filenames = sorted(glob.glob(filenames_s)) 
    lastFile = filenames[-1] 
    prefix = lastFile.split('.')[0].split('_')[0]
    timestamp = int(lastFile.split('.')[0].split('_')[-1])
    for ind in range(append): 
        timestamp += 1
        cmd = 'cp %s %s_%05u.png' %(lastFile, prefix, timestamp)
        print cmd
        call(cmd, shell=True) 

def ReadDispLine(inFile): 
    pos = np.zeros(3)
    rot = np.zeros(4)
    ts = np.fromfile(inFile, dtype=float   , count=1)
    ID = np.fromfile(inFile, dtype=np.int32, count=1)
    pos[0] = np.fromfile(inFile, dtype=float, count=1)
    pos[1] = np.fromfile(inFile, dtype=float, count=1)
    pos[2] = np.fromfile(inFile, dtype=float, count=1)
    rot[0] = np.fromfile(inFile, dtype=float, count=1)
    rot[1] = np.fromfile(inFile, dtype=float, count=1)
    rot[2] = np.fromfile(inFile, dtype=float, count=1)
    rot[3] = np.fromfile(inFile, dtype=float, count=1)
    EOL = np.fromfile(inFile, dtype=np.int32, count=1)
    # print ts, ID, pos, rot, EOL
    return pos

if __name__ == '__main__':
    ##
    demoFolder = '/home/jui-hsien/code/acoustics/work/demo/acc_noise/ball_plane'
    template = 'ball_plane_template.xml'
    tmpFilename = 'tmp.xml'
    outFile = 'out/output_%05u.png'
    cmd = 'mitsuba -DpFile=%s -o %s %s' %('%s', '%s', template)
    
    ##
    dispFile = '%s/displace.bin' %(demoFolder)
    stream = open(dispFile, 'r')
    
    count = 0
    countImg = 0
    while True:
        try:
            pos = ReadDispLine(stream)
            s = CreateObjectXML(pos[0], pos[1], pos[2])
            if (count % 5 == 0):
                outName = outFile %(countImg)
                Render(cmd, template, tmpFilename, outName, s)
                countImg += 1
            count += 1
    
        except ValueError:
            print 'Done'
            break
