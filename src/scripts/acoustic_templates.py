from material_parameters import materials

def wave_config(sim_params, objct):
    material = materials()[objct['material']]
    wave = objct['wave_params']
    return """
        <?xml version="1.0" ?>
        <impact>
        	<mesh file="%s" scale="1.0" />

          <acoustic_transfer
            distancefield="%s"
            fieldresolution="%d"
            gridresolution="%d"
            gridscale="%f"
            timestepfrequency="%d"
            substeps="%d"
            nbar="5"
            mode="0"
            radius_multipole="%f"
            mode_data_file="%s"
            rigidprefix="%s"
            rigiddensity="%f"
            outputfile="pulsedata_clipped"
            multipole_outputfile="%s"
          />

          <acceleration_pulse
            distancefield="plate"
            fieldresolution="100"
            gridresolution="300"
            OLDgridresolution="150"
            gridscale="2.95"
            OLDgridscale="0.98333"
            OLDtimestepfrequency="576000"
            timestepfrequency="440000"
            OLDsubsteps="3"
            substeps="1"
            OLDtimescale="0.000053154"
            timescale="0.000027358"
            OLDlisteningresolution="60"
            listeningresolution="40"
            listeningradius="2.0"
            numshells="5"
            numterms="2"
            radiusmultiplier="1.1892"
            rigidprefix="plateRigid"
            rigiddensity="2300.0"
            OLDoutputfile="pulsedata"
            outputfile="pulsedata_clipped"
            multitermoutputfile="pulsemapdata"
            compressedoutputfile="pulsemapdata_comp_0.04"
            compressiontolerance="0.04"
          />

        </impact>
        """ % ( objct['objFile'], objct['objFile'].replace(".obj", ""), \
                wave['field_resolution'], wave['grid_resolution'], wave['grid_scale'], \
                wave['time_step_freq'], wave['substeps'], wave['radius_multipole'], \
                objct['objPath'] + '/' + objct['objName']+".modes", objct['dataPrefix'], material['density'], \
                objct['objName']+".multipole")

def simulation_config(sim_params, objectDict, objects):
    mats = materials()

    sim_tem = """
simulation:
{
    step_size   = %f;
    time_len    = %f;
};
    """

    obj_tem = """
    {
        model = "%s";
        dx = %f;
        dy = %f;
        dz = %f;
        vx = %f;
        vy = %f;
        vz = %f;
        rot = [ %f, %f, %f, %f ];
        material = "%s";
    }"""

    mat_tem = """
    { 
        name = "%s";
        density        = %f;
        rest_coeff     = %f;
        friction_coeff = %f;
        young_modulus  = %f;
        poisson_ratio  = %f;
    }"""

    gen_tem = """
%s = (
    %s
);
    """

    finalstr = "# Generated by acoustic_templates"
    finalstr = finalstr + sim_tem%(sim_params['sim_step'], sim_params['simulation_length'])

    matl = []
    for material in mats:
        m = mats[material]
        matl.append(mat_tem%(material, m['density'], m['restitution_coeff'], m['friction_coeff'],\
                            m['youngModulus'], m['poissonRatio']))

    finalstr = finalstr + gen_tem%('materials', ',\n'.join(matl))

    objl = []
    for obj in objects:
        p = obj['position']
        v = obj['velocity']
        r = obj['rotation']
        objl.append(obj_tem%(objectDict[obj['object']]['tetFile'], \
                             p[0], p[1], p[2], \
                             v[0], v[1], v[2], \
                             r[0], r[1], r[2], r[3], \
                             objectDict[obj['object']]['material']))
    finalstr = finalstr + gen_tem%('objs', ',\n'.join(objl))

    return finalstr

def audio_gen_config(sim_params, objectDict, objects):
    temp = """
# Generated by acoustic_templates
sndrate   = %d;
sim_step  = %f;
sound_len = %f;

# the file storing impulses
impulse_file = "%s";
scale_prefix = "%s";

sndobj = (
    %s
);
    """
    obj_tem = """
    { 
        geometry       = "%s";
        eigenmodes     = "%s";
        transfer_file  = "%s";
        tick_file      = "%s";
        alpha          = %.12f;
        beta           = %f;
        id             = %d;
        dampingGamma   = %.12f ;
        ear_angle      = "ear_angles-%d.txt";
        numfixed       = 0; 
        density        = %f;
        contact_damping_winlen    = 800;
        contact_damping_threshold = 800;
    }"""

    objl = []
    for i in range(0, len(objects)):
        obj = objects[i]
        typ = objectDict[obj['object']]
        m = materials()[typ['material']]

        objl.append(obj_tem%(typ['objPath']+'/'+typ['objName']+'.geo.txt', typ['objPath']+'/'+typ['objName']+'.modes',\
                             sim_params['transferFilePattern'].replace('%d', str(i), 1),\
                             sim_params['outputFilePattern']%(i),\
                             m['beta'], m['alpha'], i, m['gamma'], i, m['density']
                             ))

    return temp%(sim_params['sound_rate'], sim_params['sim_step'], sim_params['simulation_length'],\
                 sim_params['modalImpulse'], sim_params['forceScalePrefix'], ",\n".join(objl))

if __name__ == '__main__':
    import params
    print audio_gen_config(params.sim_parameters(), params.objectDict(), params.objects())

def estimate_timescales_gen(sim_params, objectDict, objects):
    template = """<?xml version="1.0" ?>
<impact>
  <scene>
    <meshList>
      %s
    </meshList>

    <objectList>
      %s
    </objectList>
  </scene>
</impact>"""

    mesh_tem = """<mesh
        file="%s"
        id="%s"
        rigidprefix="%s"
        distancefield="%s"
        fieldresolution="%d"
        OLDpulsefile="../coin_mesh/old_pulsedata/pulsedata"
        pulsefile="../plate/pulsemapdata"
        density="%f"
      />"""

    object_tem = """<object id="%s" />"""

    meshes = []
    for objType in objectDict:
        typ = objectDict[objType]
        m = materials()[typ['material']]
        meshes.append(mesh_tem%(typ['objFile'], typ['objName'], typ['dataPrefix'], \
                                typ['objFile'].replace(".obj", ""), \
                                typ['wave_params']['field_resolution'], m['density']))

    objects_t = []
    for obj in objects:
        objects_t.append(object_tem%(obj['object']))

    return template%("\n".join(meshes), "\n".join(objects_t))