#!/usr/bin/env python
import sys
from get_paths import *

def SplitObjects(fileLines): 
    objectLinePos = []
    objectLines = []
    objectNames = []
    objectNumVertices = []
    lineNum = 0
    for line in fileLines: 
        strings = line.split()
        if strings[0] == 'o': 
            objectLinePos.append(lineNum)
            objectNames.append(strings[1]) 
        lineNum += 1
    N_objects = len(objectLinePos)
    count = 0
    accumCountVert = 0
    for obj_lp in objectLinePos: 
        if obj_lp != objectLinePos[-1]:
            lineRange = [obj_lp, objectLinePos[count+1]]
        else: 
            lineRange = [obj_lp, len(fileLines)]
        objectLines.append(fileLines[lineRange[0]:lineRange[1]])
        
        # count num of vertices and renumber face vertices
        N_lines = lineRange[1] - lineRange[0]
        count_v = 0
        for l in objectLines[len(objectLines)-1]:
            strings = l.split()
            if strings[0] == 'v':
                count_v += 1
          
        count_l = 0
        for l in objectLines[len(objectLines)-1]: 
            strings = l.split()
            if strings[0] == 'f': 
                N_faceVert = len(strings) - 1 
                renumberedVert = 'f '
                for fv_str in strings[1::]: 
                    fv_num = int(fv_str) - accumCountVert # now everything should start from 0
                    renumberedVert += '%u ' %(fv_num)
                renumberedVert += '\n'

                objectLines[len(objectLines)-1][count_l] = renumberedVert
            count_l += 1
        
        accumCountVert += count_v
        count += 1

    print '%u objects read.' %(len(objectNames))
    return objectNames, objectLines

def WriteObjectsToFiles(objectNames, objectLines): 
    N_objects = len(objectNames)
    assert(N_objects == len(objectLines))

    objFilenames = []
    for o_idx in range(N_objects):
        objName = objectNames[o_idx]
        objLines = objectLines[o_idx]
        filename = '%s.obj' %(objName)
        print 'Writing object %s to file %s' %(objName, filename)
        of = open(filename, 'w')
        for line in objLines:
            of.write(line)
        of.close()
        objFilenames.append(filename)
    return objFilenames

def PrintAsConfigFormat(absDirPath, objFilenames):
    print 'START print config segment'
    count = 0
    for obj in objFilenames:
        objNameNoAffix = obj[:-4]
        rigid_object_str = """
           <rigid_object
               working_directory="%s"
               object_prefix="%s"
               id="%u"
               fieldresolution="400"
               scale="1.0"
               initial_position_x="0.00000000"
               initial_position_y="0.00000000"
               initial_position_z="0.00000000"
               material_id="0"
               /> """ %(absDirPath, objNameNoAffix, count)
        print rigid_object_str
        count += 1
    print 'END print config segment'

if __name__ == '__main__': 
    if len(sys.argv) < 2: 
        print '**Usage: %s <filename> [absPath]' %(sys.argv[0])
        sys.exit()

    if len(sys.argv) == 3:
        absPath = sys.argv[2]
    else: 
        absPath = projectPath() + '/work/pool_cooler'

    filename = sys.argv[1]
    fileLines = open(filename).readlines()
    objNames, objLines = SplitObjects(fileLines)
    objFilenames = WriteObjectsToFiles(objNames, objLines)
    PrintAsConfigFormat(absPath, objFilenames)

