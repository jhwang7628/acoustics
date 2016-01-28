#!/usr/bin/python
#
# Script which runs the initial rigid-body simulation

import sys
import os
import config_path

def encode_frames(sim_params):
	encodeScript = '%s/src/scripts/encode_png_sequence.py' % ( config_path.installPath() );

	# We ran a rigid-body simulation at 1000Hz, then dumped out every 10th time
	# step, resulting in a rate of 100 FPS
	frameRate = 1/(sim_params['sim_step']*sim_params['frameFrequency']);

	imageXres = 1024;
	imageYres = 576;

	# Make a movie!
	cmd = 'python %s %s/frame no_audio_%s %d %d %d' \
	       % ( encodeScript, sim_params['frameDir'], sim_params['videoFile'], imageXres, imageYres, frameRate );

	os.system(cmd);
