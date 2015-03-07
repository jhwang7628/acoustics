#!/usr/bin/python
#
# Script which runs the initial rigid-body simulation

import sys
import os
import config_path
import acoustic_templates

def run_simulation(sim_params, objectDict, objects):
	with open(sim_params['simulationFile'], 'w+') as sFile:
		sFile.write(acoustic_templates.simulation_config(sim_params, objectDict, objects))

	cmd = '%s/rigidsim %s' % ( config_path.binPath(), sim_params['simulationFile']);
	print( "=== Running the rigid-body simulator" );
	os.system(cmd);
	print( "=== Done" );
