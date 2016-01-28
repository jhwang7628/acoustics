#!/usr/bin/python
#
# Script which runs the initial rigid-body simulation

import sys
import os
import config_path

def make_frames(sim_params):
	# Put the frames in ./frames
	os.system('mkdir -p %s' % ( sim_params['frameDir'] ) );

	# Bring up the replayer
	cmd = '%s/replayer -f %s -d %s -r %d -o %s/frame' \
	       % ( config_path.binPath(), sim_params['simulationFile'], sim_params['displaceFile'], sim_params['frameFrequency'], sim_params['frameDir']);

	print( "=== Runing the rigid-body replayer" );
	print( "=== To generate frames:" );
	print( "===     1) Adjust the camera view as desired" );
	print( "===     2) Press 'P' to enable frame snapshots" );
	print( "===     3) Press space to start writing frames" );
	print( "" );

	os.system(cmd);
