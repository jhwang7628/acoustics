#!/usr/bin/env python
import struct
from scene import *
from get_paths import *
################################################################################
## Format sees Simulation.cpp:Write_Output_Alltransform
################################################################################
def Read_Transform_Frames(filename, rigid_objects): 
    for o in rigid_objects: 
        assert(o.restCOM is not None)
    ifs = open(filename,'rb')
    all_bytes = ifs.read()
    p = 0
    count_frames = 0
    while p < len(all_bytes):
        frame_id = struct.unpack('d', all_bytes[p:p+8])[0]
        p += 8
        count_frames += 1
        count_objects_transforms = 0
        while True: 
            obj_id = struct.unpack('i', all_bytes[p:p+4])[0]
            p += 4
            if obj_id == -1:
                break
            x0 = rigid_objects[obj_id].restCOM
            x  = struct.unpack('ddd', all_bytes[p:p+24])
            p += 24
            quaternion = struct.unpack('dddd', all_bytes[p:p+32])
            rotAxis, rotDeg = Rotation.Quaternion_To_Axis_Rotation_Degree(quaternion)
            p += 32
            frame = Rigid_Frame(frame_id)
            frame.Add_Translation([-x0[0], -x0[1], -x0[2]])
            frame.Add_Rotation(rotAxis, rotDeg)
            frame.Add_Translation(x)
            rigid_objects[obj_id].frames.append(frame)
            count_objects_transforms += 1
        assert(count_objects_transforms == len(rigid_objects))
    ifs.close()
    print 'Frames read from file %s: %d' %(filename, count_frames)
    return rigid_objects

def Read_Object_Kin_Frames(filename, rigid_object, fix_blender_rotation=True):
    '''
    Read in tranformation info for all frames from .kin file for a single rigid object.
    '''
    with open(filename, 'r') as fp:
        for line in fp:
            tokens = line.strip().split()
            if len(tokens) > 0:
                frame_id = int(tokens[0])
                translation = [float(tokens[1]), float(tokens[2]), float(tokens[3])]
                quat = [float(tokens[4]), float(tokens[5]), float(tokens[6]), float(tokens[7])]
                # clamp w to avoid numerical error
                quat[0] = max(-1.0, min(1.0, quat[0]))
                rotAxis, rotDeg = Rotation.Quaternion_To_Axis_Rotation_Degree(quat)

                print 'Frame: {}, Trans: {}, Quat: {}, RotAxis: {}, Deg: {}'.format(frame_id, translation, quat, rotAxis, rotDeg)

                # transforms to object for this frames
                frame = Rigid_Frame(frame_id)
                if fix_blender_rotation:
                    # blender uses z-up as default so we apply -90 rot about x-axis
                    frame.Add_Rotation([1.0, 0.0, 0.0], -90.0)
                frame.Add_Rotation(rotAxis, rotDeg)
                frame.Add_Translation(translation)
                rigid_object.frames.append(frame)

    return rigid_object

def Read_Object_Kin_Frames_Upsample(filename, rigid_object, input_fps, new_fps, total_frames, fix_blender_rotation=True):
    '''
    Read in tranformation info for all frames from .kin file for a single rigid object.
    Upsamples the frame data in the .kin file from init_fps to new_fps to use for
    an animation (w/ the new_fps) that is total_frames in length to be rendered. Uses linear interp
    for translation and slerp for rotation.
    '''
    input_interval = 1.0 / float(input_fps)
    with open(filename, 'r') as fp:
        # get first 2 lines for initial frames
        left_frame_id, left_trans, left_quat = Read_Next_Line_Kin(fp)
        right_frame_id, right_trans, right_quat = Read_Next_Line_Kin(fp)
        for cur_frame in range(0, total_frames):
            # check if passed current right frame
            cur_time = float(cur_frame) / float(new_fps)
            diff = ((float(right_frame_id) - 1.0) / float(input_fps)) - cur_time
            alpha = 1.0 - diff / input_interval
            if diff < 0:
                # update adjacent data
                left_frame_id, left_trans, left_quat = right_frame_id, list(right_trans), list(right_quat)
                right_frame_id, right_trans, right_quat = Read_Next_Line_Kin(fp)
                if right_frame_id == -1:
                    right_frame_id = left_frame_id
                    right_trans = left_trans
                    right_quat = left_quat
                # set alpha with updated data
                diff = ((float(right_frame_id) - 1.0) / float(input_fps)) - cur_time
                alpha = 1.0 - diff / input_interval
            # interpolate
            translation = Translation.Lerp(left_trans, right_trans, alpha)
            quat = Rotation.Quaternion_Slerp(left_quat, right_quat, alpha)
            rotAxis, rotDeg = Rotation.Quaternion_To_Axis_Rotation_Degree(quat)

            print 'Frame: {}, Trans: {}, Quat: {}, RotAxis: {}, Deg: {}'.format(cur_frame+1, translation, quat, rotAxis, rotDeg)

            # transforms to object for this frames
            frame = Rigid_Frame(cur_frame+1)
            if fix_blender_rotation:
                # blender uses z-up as default so we apply -90 rot about x-axis
                frame.Add_Rotation([1.0, 0.0, 0.0], -90.0)
            frame.Add_Rotation(rotAxis, rotDeg)
            frame.Add_Translation(translation)
            rigid_object.frames.append(frame)
    return rigid_object

def Read_Next_Line_Kin(fp):
    ''' Reads line from a kin file and returns frame_id, translation vect, and quaternion. '''
    frame_id = -1
    translation = []
    quat = []
    line = fp.readline()
    tokens = line.strip().split()
    if len(tokens) > 0:
        frame_id = int(tokens[0])
        translation = [float(tokens[1]), float(tokens[2]), float(tokens[3])]
        quat = [float(tokens[4]), float(tokens[5]), float(tokens[6]), float(tokens[7])]
        # clamp w to avoid numerical error
        quat[0] = max(-1.0, min(1.0, quat[0]))
    return frame_id, translation, quat


if __name__ == "__main__":
    ##
    test_file = projectPath() + '/work/demo/wobbling_plate/displace.bin'
    plate = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/plate/plate.obj')
    plate.Read_Rest_COM('/home/jui-hsien/data/models/plate/plate_centerOfMass.3vector')
    ##
    rigid_objects = [plate]
    Read_Transform_Frames(test_file, rigid_objects)

