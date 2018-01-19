#!/usr/bin/env python
import struct,glob,re
import numpy as np
from get_paths import *
################################################################################
################################################################################
class Shell_Frame:
    def __init__(self, owner, ID, disp):
        self.owner = owner
        self.ID = ID
        self.disp = disp
    def Format_String(self): # TODO
        return ''

################################################################################
################################################################################
def Read_Shell_Disp(directory, shellobj, ext='.displacement',
                    read_constraints=False):
    print 'Read_Shell_Disp from dir: %s' %(directory)
    if shellobj.match_re is None:
        filenames = sorted(glob.glob(directory+'/*%s' %(ext)))
    else:
        filenames = [directory+'/'+f for f in os.listdir(directory) if re.search(shellobj.match_re, f)]
        filenames = sorted(filenames)

    for f in filenames:
        print '  Reading shell displacement: ', f
        idx = int(f.split('/')[-1].split('.')[0])
        with open(f, 'rb') as stream:
            N = struct.unpack('i', stream.read(4))[0]
            disp = struct.unpack('d'*N, stream.read(8*N))
            frame = Shell_Frame(shellobj, idx, disp)
            shellobj.frames[idx] = frame
    if read_constraints:
        # TODO
        print 'read constraints'

################################################################################
################################################################################
def Read_Shell_Vertex_Map(directory, shellobj, filename='vertex_map.txt'):
    print 'Read_Shell_Vertex_Map from dir: %s' %(directory)
    f = glob.glob('%s/%s' %(directory, filename))[0]
    shellobj.v_map = np.loadtxt(f, dtype=int)

################################################################################
################################################################################
def Read_Shell_Render2Sim_Bary_Map(shellobj, filename):
    print 'Read_Shell_Render2Sim_Bary_Map: ', filename
    with open(filename, 'r') as stream:
        lines = stream.readlines()
        for idx in range(len(lines)):
            if lines[idx][0] != '#':
                break
        for ii in range(idx, len(lines)):
            tokens = lines[ii].split()
            tri_id = int(tokens[0])
            bary = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])
            shellobj.render2sim_bary_map.append([tri_id, bary])

if __name__ == "__main__":
    ##
    test_file = projectPath() + '/work/demo/wobbling_plate/displace.bin'
    plate = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/plate/plate.obj')
    plate.Read_Rest_COM('/home/jui-hsien/data/models/plate/plate_centerOfMass.3vector')
    ##
    rigid_objects = [plate]
    Read_Transform_Frames(test_file, rigid_objects)

