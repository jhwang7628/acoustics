#!/usr/bin/python
#
# Estimates impact timescales using Hertz contact theory

import sys
import os
import config_path
import acoustic_templates

def estimate_timescales(sim_params, objectDict, objects):

	impulseFile = sim_params['impulseFile'];
	forceScalePrefix = sim_params['forceScalePrefix'];
	with open("default.xml", 'w+') as sFile:
		sFile.write(acoustic_templates.estimate_timescales_gen(sim_params, objectDict, objects))

	cmd = '%s/estimate-timescales %s %s' \
	       % ( config_path.binPath(), impulseFile, forceScalePrefix );

	print( "=== Estimating Hertz contact time scales" );
	os.system(cmd);
	print( "=== Done" );
