#!/usr/bin/env python
import sys,os
import numpy as np
sys.path.insert(0, '/home/jui-hsien/code/acoustics/src/scripts/rendering')
from rigid_body_data import Read_Transform_Frames,Read_Object_Kin_Frames_Upsample
from scene import Object,Rigid_Wavefront_Obj,Ply,Shell_Obj
from scene_templates import Simple_3PointLight_Diffuse_Mts_Scene
from subprocess import call
################################################################################
################################################################################
if __name__ == "__main__":
    ## Input
    input_kin_fps = 1290
    output_fps = 44100
    output_directory = 'frames'
    total_time = 2.165
    shell_data_front_pad_frames = 208
    shell_data_front_pad_time = float(shell_data_front_pad_frames)/float(input_kin_fps)
    # frames = range(0, 88200+int(shell_data_front_pad_time*output_fps), 735)
    frames = [7144, 8085, 12495]

    obj_cymbal = Shell_Obj('../assets/cymbal/cymbal_aligned_assets_18inch_cymbal-only.obj')
    # obj_stick  = Rigid_Wavefront_Obj('../assets/drumsticks/drumsticks_scaled.obj', use_COM=False)
    # obj_stand  = Rigid_Wavefront_Obj('../assets/cymbal/cymbal_aligned_assets_18inch_stand-only.obj', use_COM=False)
    shell_objects = [obj_cymbal]
    # dynamic_objects = [obj_stick]
    # static_objects = [obj_stand]
    dynamic_objects = []
    static_objects = []
    # kin_files = ['../assets/drumstick_2_y-up.kin']
    shell_data_dir = '../imp_fullspace_sim-011'
    render_sim_map = 'test.map'

    ## Manual settings
    mts_bin_template = 'mitsuba %s -o %s'
    integrator_string = """
    <integrator type="photonmapper">
        <integer name="directSamples" value="16"/>
        <integer name="glossySamples" value="32"/>
    </integrator>
    """
    obj_cymbal.Set_Material_String(None)
    obj_stick.Set_Material_String(None)
    obj_stand.Set_Material_String("""
    <bsdf name="Metal.001" type="roughconductor">
        <string name="material" value="Al"/>
    </bsdf>
    <bsdf name="Black_Plastic.001" type="roughplastic">
        <rgb name="diffuseReflectance" value="0.02, 0.02, 0.02"/>
    </bsdf>""")

    ## Derived
    total_frames = int(total_time*output_fps)

    ## Start
    if not os.path.isdir(output_directory):
        call('mkdir -p %s' %(output_directory), shell=True)
    # for i, obj in enumerate(dynamic_objects):
    #     Read_Object_Kin_Frames_Upsample(kin_files[i], obj, input_kin_fps, output_fps, total_frames)
    count = 0
    for frame_id in frames:
        print '\n\nRendering frame %u\n' %(frame_id)
        shell_frame_id = max(0, frame_id - int(shell_data_front_pad_time*output_fps))

        frame_str = r'%0.9u.displacement' %(shell_frame_id)
        obj_cymbal.Set_MatchRegex(frame_str)
        obj_cymbal.Read_Data(shell_data_dir, render_sim_map)

        frame = obj_cymbal.frames[shell_frame_id]
        # frame = obj_cymbal.frames[key]
        output_file = '%s/%s-%.5d.png' %(output_directory, 'cymbal_heavy', count)
        if os.path.isfile(output_file):
            count += 1
            continue
        frame_string = ''
        for o in shell_objects:
            frame_string += o.Format_String(frame)
        for o in dynamic_objects:
            frame_string += o.Format_String(frame_id)
        for o in static_objects:
            frame_string += o.Format_String(-1)

        camera_origin = [-1.80351, 0.752225, 0.831354]
        # close up
        fov = 18
        camera_target = [0.,0.716, 0. ]
        # overview
        # fov = 60
        # camera_target = [0.,0.60, 0. ]
        config_string = Simple_3PointLight_Diffuse_Mts_Scene(frame_string,
                                                             definition='test',
                                                             fov=fov,
                                                             key_light_value=80,
                                                             fill_light_value=70,
                                                             rim_light_value=15,
                                                             camera_x=camera_origin[0],
                                                             camera_y=camera_origin[1],
                                                             camera_z=camera_origin[2],
                                                             lookat_x=camera_target[0],
                                                             lookat_y=camera_target[1],
                                                             lookat_z=camera_target[2],
                                                             # integrator_string=integrator_string,
                                                             sample_count=32,
                                                             ground_rotation_deg=-90.,
                                                             fill_light_origin=[2.,1.,7],
                                                             rim_light_origin=[6.,2.,-2.],
                                                             scene_height=0.7
                                                             )

        # write strings to config
        tmp_xml = 'tmp_%s.xml' %('test')
        ofs = open(tmp_xml, 'w')
        for s in config_string:
            ofs.write(s)
        ofs.close()
        cmd = mts_bin_template %(tmp_xml, output_file)
        print cmd
        call(cmd, shell=True)
        count += 1
