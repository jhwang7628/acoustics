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
            timestepfrequency="192000"
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
                objct['objName']+".modes", objct['dataPrefix'], material['density'], \
                objct['objName']+".multipole")


