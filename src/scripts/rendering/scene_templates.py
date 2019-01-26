#!/usr/bin/env python
################################################################################
################################################################################
def Simple_3PointLight_Diffuse_Mts_Scene(frame_string, definition='test', resolution=None, integrator_string=None, fov=45.0, ground_offset_y=0.0, camera_x=0.05, camera_y=1.0, camera_z=4.0, lookat_x=0.0, lookat_y=0.5, lookat_z=0.0, sample_count=60, key_light_value=50, fill_light_value=18, rim_light_value=15, ground_material='default',camera_up=[0,1,0], ground_rotation_axis=[0,1,0], ground_rotation_deg=0.0, fill_light_origin=[-2, 1, 7], rim_light_origin=[0, 3, -9], scene_height=0.0):
    ####################################
    ## default settings
    if resolution is None:
        aspect_ratio=16.0/9.0
        definitions={'480p' : [480. ,int(480./aspect_ratio)],
                     'test' : [960. ,int(960./aspect_ratio)],
                     '4k'   : [3840.,int(3840./aspect_ratio)],
                     '1080p': [1080.,int(1080./aspect_ratio)],
                     'HD'   : [1920.,int(1920./aspect_ratio)]}
        resolution = definitions[definition]
    if ground_material == 'diffuse_white':
        ground_material_id = 'diffuse_white_ground'
    else: 
        ground_material_id = 'default_ground'

    if integrator_string is None:
        integrator_string = """
    <integrator type="path">
        <integer name="maxDepth" value="120"/>
    </integrator> """

    ####################################
    s = """<scene version="0.5.0">
    %s
    <bsdf type="diffuse" id="default_ground">
        <spectrum name="reflectance" value="0.2 0.2 0.2"/>
    </bsdf>
    <bsdf type="diffuse" id="diffuse_material">
        <spectrum name="reflectance" value="0.2 0.2 0.2"/>
    </bsdf>
    <bsdf type="diffuse" id="diffuse_white_ground">
        <spectrum name="reflectance" value="1.0 1.0 1.0"/>
    </bsdf> 
    %s
    <shape type="obj">
        <string name="filename" value="/home/jui-hsien/data/models/render_assets/infinite_ground.obj"/>
        <transform name="toWorld">
            <rotate x="%f" y="%f" z="%f" angle="%f"/>
            <translate x="0.0" y="%f" z="0.0"/>
        </transform>
        <ref id="%s"/>
    </shape>
    <bsdf type="diffuse" id="light">
        <spectrum name="intensity" value="500"/>
    </bsdf>
    <!-- key light -->
	<shape type="obj">
        <string name="filename" value="render_assets/area_light.obj"/>
        <transform name="toWorld">
            <scale value="1.0"/>
            <!-- this rotation takes vector (0,-1,0) to (-1,-1,-1)-->
            <rotate x="0.70710678" y="0." z="-0.70710678" angle="54.735610317245346"/>
            <translate x="3.0" y="3.0" z="3.0"/>
            <rotate x="%f" y="%f" z="%f" angle="%f"/>
        </transform>
        <ref id="light"/>
        <emitter type="area">
            <spectrum name="radiance" value="%f"/>
        </emitter>
    </shape>
    <!-- fill lights -->
    <emitter type="spot">
        <transform name="toWorld">
            <lookat origin="%f %f %f" target="0 %f 0"/>
        </transform>
        <float name="cutoffAngle" value="40"/>
        <spectrum name="intensity" value="%f"/>
    </emitter>
    <!-- rim light -->
    <emitter type="spot">
        <transform name="toWorld">
            <lookat origin="%f %f %f" target="0 %f 0"/>
        </transform>
        <float name="cutoffAngle" value="30"/>
        <spectrum name="intensity" value="%f"/>
    </emitter>
    <sensor type="perspective">
        <transform name="toWorld">
            <lookat origin="%f %f %f"
                    target="%f %f %f"
                    up="%f, %f, %f"/>
        </transform>
        <float name="fov" value="%f"/>
        <string name="fovAxis" value="x"/>
        <film type="ldrfilm">
            <integer name="width" value="%d"/>
            <integer name="height" value="%d"/>
            <string name="fileFormat" value="png"/>
            <string name="pixelFormat" value="rgba"/>
            <boolean name="banner" value="false"/>
        </film>
        <sampler type="ldsampler">
            <integer name="sampleCount" value="%d"/>
        </sampler>

    </sensor>
</scene> """ %(integrator_string, frame_string, ground_rotation_axis[0], ground_rotation_axis[1], ground_rotation_axis[2], ground_rotation_deg, ground_offset_y, ground_material_id, ground_rotation_axis[0], ground_rotation_axis[1], ground_rotation_axis[2], ground_rotation_deg, key_light_value, fill_light_origin[0], fill_light_origin[1], fill_light_origin[2], scene_height, fill_light_value, rim_light_origin[0], rim_light_origin[1], rim_light_origin[2], scene_height, rim_light_value, camera_x, camera_y, camera_z, lookat_x, lookat_y, lookat_z, camera_up[0], camera_up[1], camera_up[2], fov, resolution[0], resolution[1], sample_count)
    return s

