#!/usr/bin/python
#
# Script which runs the initial rigid-body simulation

import sys
import os
import config_path

def prep_transfer_eval(sim_params, objectDict, objects):

	nObjects = len(objects);

	listenPosition = sim_params['listeningPosition']

	# Evaluate transfer at 1000Hz
	transferEvalTimestep = 1.0/sim_params['transferFrequency']

	# Build mc.txt
	print( "=== Preparing mass center file" );
	massCenterDict = {}
	tempfile = "_temp_mc_.txt"
	for objectType in objectDict:
		os.system( 'rm -f %s' % tempfile );
		cmd = '%s/init-mass-centers %s 0 0 %s' % (	config_path().binPath() \
													objectDict[objectType]['tetFile'],\
													tempfile );
		os.system(cmd);
		with open(tempfile) as tfile:
			for line in tfile:
				w = line.split()
				x = float(w[1])
				y = float(w[2])
				z = float(w[3])
				massCenterDict[objectType] = (x, y, z)

	os.system( 'rm -f %s' % tempfile );
	os.system( 'rm -f %s' % sim_params['massCenterFile']);
	with open(sim_params['massCenterFile'], 'w+')	as outFile:
		for i in range(0, len(objects)):
			obje = objects[i]
			cm = massCenterDict[obje['object']]
			outFile.write("%d %f %f %f\n"%(i, cm[0], cm[1], cm[2]))

	print( "=== Done\n" );

	print( "=== Preparing object <--> model map file" );

	# Build obj_model_map.txt (used by transfer scripts on the hydras)
	os.system( 'rm -rf %s' % ( sim_params['mapFile'] ) );

	with open(sim_params['mapFile']) as mapFile:
		for i in range(0, len(objects)):
			obje = objects[i]
			mapFile.write("%d %s\n"%(i, obje['object']))
	print( "=== Done" );

	# Generate the fieldpoints file(s)

	print( "=== Setting up field points" );
	os.system( 'mkdir -p field_pts' );
	cmd = '%s/fieldpoints_gen -i %s -x %s -o %s -n %d -t %f %f %f %f' \
	       % ( 	config_path().binPath(), sim_params['displaceFile'], \
	       		sim_params['massCenterFile'], sim_params['outputFilePattern'], \
	          	nObjects, transferEvalTimestep, \
	           	listenPosition[0], listenPosition[1], listenPosition[2] );

	os.system(cmd);
	print( "=== Done" );

	# simName = 'plate_drop_test';

	# # Finally, put everything in a tarball so that it can be moved to the cluster
	# # 
	# # Move plate_drop_test.tar.gz t= /hydra/S1/impact_snd/fastbem/sim and untar
	# # it there.  This should create the directory 'plate_drop_test'
	# #
	# # Navigate to /hydra/S1/impact_snd/fastbem/sim/plate_drop_test
	# # From this directory, run ./../../scripts/comp_trans_val.pl . <objid> <modeid>
	# # For this example, objid will be 0, since there is only one object, and
	# # you should run this for all audible modes.
	# #
	# # Once you are done with all of the modes, copy
	# # /hydra/S1/impact_snd/fastbem/sim/plate_drop_test/fbem_ret
	# # back to the plate_drop_test work directory on your machine.  This directory
	# # will store files with time series of values for each modal transfer function
	# cmd = "tar --transform 's,^,%s/,' -czvf %s.tar.gz field_pts/ %s" \
	#        % ( simName, simName, mapFile );
	# os.system(cmd);
