#!/usr/bin/env python
import sys,os
import numpy as np
sys.path.insert(0, '/home/antequ/code/acoustics/src/scripts/rendering')
from rigid_body_data import Read_Transform_Frames
from scene import Rigid_Wavefront_Obj,Ply
from scene_templates import Simple_3PointLight_Diffuse_Mts_Scene
from subprocess import call
from get_paths import *
################################################################################
################################################################################
if __name__ == "__main__":
    ##
    rbd_file = projectPath() + '/work/demo/acc_noise/ball_bowl_ceramics/displace.bin'
    ##
    output_directory = 'frames'
    mts_bin_template = 'nice -n 2 mitsuba %s -o %s'
    integrator_string = """
    <integrator type="photonmapper">
        <integer name="directSamples" value="16"/>
        <integer name="glossySamples" value="32"/>
    </integrator>
    """
    ##
    if not os.path.isdir(output_directory):
        call('mkdir -p %s' %(output_directory), shell=True)
    obj_ball = Rigid_Wavefront_Obj('/home/jui-hsien/data/models/balls/ball_d2cm.obj')
    obj_mug  = Rigid_Wavefront_Obj('/home/antequ/data/models/bowl_clean_ceramics/bowl.obj')
    obj_ball.Set_Material_String("""
    <bsdf type="roughconductor" id="sphere">
        <string name="material" value="Al"/>
    </bsdf>""")
    obj_mug.Set_Material_String("""
    <bsdf type="dielectric" id="glass">
        <string name="intIOR" value="bk7"/>
        <string name="extIOR" value="air"/>
        <float name="maxSmoothAngle" value="0"/>
    </bsdf>""") 
    rigid_objects = [obj_ball, obj_mug]
    Read_Transform_Frames(rbd_file, rigid_objects)
    framesstart = 0.0
    framesend = 10000.0
    total_time = 1.0
    fps = 60
    # frames = [2000]
    # frames = range(0, 10000, 10000.0/30.0)
    frames = np.around(np.arange(framesstart,framesend+0.5,framesend/(total_time * fps))).astype(int)
   # frames = [10000]
    for frame in frames: 
        output_file = '%s/%s-%.5d.png' %(output_directory, 'test', frame)
        if os.path.isfile(output_file):
            continue
        frame_string = ''
        for o in rigid_objects:
            frame_string += o.Format_String(frame)

        camera_origin = [0.,0.3 , 4.5]
        camera_target = [0.,0.1, 0. ]
        config_string = Simple_3PointLight_Diffuse_Mts_Scene(frame_string, definition='HD', fov=10, key_light_value=80,
                                                             camera_x=camera_origin[0], camera_y=camera_origin[1], camera_z=camera_origin[2], 
                                                             lookat_x=camera_target[0], lookat_y=camera_target[1], lookat_z=camera_target[2], 
                                                             integrator_string=integrator_string,
                                                             sample_count=70
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



