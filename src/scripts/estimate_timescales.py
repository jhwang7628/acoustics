#!/usr/bin/python
#
# Estimates impact timescales using Hertz contact theory

import sys
import os
import config_path

def estimate_timescales(sim_params):

	impulseFile = sim_params['impulseFile'];
	forceScalePrefix = sim_params['forceScalePrefix'];

	cmd = '%s/estimate-timescales %s %s' \
	       % ( config_path.binPath(), impulseFile, forceScalePrefix );

	print( "=== Estimating Hertz contact time scales" );
	os.system(cmd);
	print( "=== Done" );
