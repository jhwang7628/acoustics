#!/usr/bin/python
#
# Generates sound from rigid-body simulation data

import sys
import os
import config_path
import shutil

def generate_sound(sim_params, objectDict, objects):
	with open(sim_params['soundConfigFile'], 'w+') as sFile:
		sFile.write(acoustic_templates.audio_gen_config(sim_params, objectDict, objects))

	cmd = '%s/sound-generator %s' % (config_path.binPath(), sim_params['soundConfigFile']);

	print( "=== Running modal sound generator" );
	os.system(cmd);
	print( "=== Done" );

	# Use octave to turn the output in to a wave file
	cmd = "octave --eval \"addpath('%s'); make_sound();\"" % ( config_path.matlabPath() );
	os.system(cmd);

	shutil.copy2('test.wav', sim_params['audioFile'])
