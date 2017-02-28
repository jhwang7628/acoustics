#!/usr/bin/env python
import struct
from scene import *
################################################################################
## Format sees Simulation.cpp:Write_Output_Alltransform
################################################################################
def Read_Transform_Frames(filename, rigid_objects): 
    ifs = open(filename,'rb')
    all_bytes = ifs.read()
    p = 0
    count_frames = 0
    while p < len(all_bytes):
        frame_id = struct.unpack('d', all_bytes[p:p+8])[0]
        p += 4
        count_frames += 1
        count_objects_transforms = 0
        while True: 
            obj_id = struct.unpack('i', all_bytes[p:p+4])[0]
            p += 4
            if obj_id == -1:
                break
            x0 = struct.unpack('ddd', all_bytes[p:p+24])
            p += 24
            x  = struct.unpack('ddd', all_bytes[p:p+24])
            p += 24
            rotation  = struct.unpack('dddd', all_bytes[p:p+32])
            p += 32
            frame = Rigid_Frame(frame_id)
            frame.Add_Translation([-x0[0], -x0[1], -x0[2]])
            frame.Add_Rotation(rotation[:3], rotation[3])
            frame.Add_Translation(x)
            rigid_objects[obj_id].frames.append(frame)
            count_objects_transforms += 1
        assert(count_objects_transforms == len(rigid_objects))
    ifs.close()
    print 'Frames read from file %s: %d' %(filename, count_frames)
    return rigid_objects

if __name__ == "__main__":
    ##
    test_file = '/home/jui-hsien/data/bounce_map/BOUNCE/fix/common/all-objects.transform'
    B = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/B/B.obj')
    O = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/O/O.obj')
    U = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/U/U.obj')
    N = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/N/N.obj')
    C = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/C/C.obj')
    E = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/E/E.obj')
    ##
    rigid_objects = [B, O, U, N, C, E]
    print rigid_objects
    Read_Transform_Frames(test_file, rigid_objects)

