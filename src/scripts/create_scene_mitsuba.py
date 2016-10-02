#!/usr/bin/env python
import sys
import glob

def GetHeader(): 
    return """<?xml version="1.0" encoding="utf-8"?>
<scene version="0.5.0"> """

def GetFooter(): 
    return """
</scene> """

def GetObjectToScene(objName):
    shapeTemplate = """
    <shape type="obj">
        <string name="filename" value="%s"/>
    </shape> """ %(objName)

    return shapeTemplate

if __name__ == '__main__': 
    filenames = sorted(glob.glob('*.obj')) 
    of = open('test_render_objects.xml', 'w')
    strBuffer = []
    strBuffer += GetHeader()
    for f in filenames: 
        strBuffer += GetObjectToScene(f)
    strBuffer += GetFooter()
    for s in strBuffer:
        of.write(s)
    of.close()



