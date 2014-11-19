#!/usr/bin/python
#
# Runs (some of) the steps needed to precompute modal sound data for an object

import sys
import os
from material_parameters import materials
import config_path

def precompute_modes(object_param):

	# By default, don't product output, just write the commands
	redirect = '/dev/null';
	#redirect = '/dev/stdout';

	installPath = config_path.installPath()
	binPath = config_path.binPath()
	matlabPath = config_path.matlabPath()

	material = materials()[object_param['material']]

	print("=== Running the isostuffer ===")
	cmd = '%s/isostuffer -R %d -L %d -M %d -a %f -b %f %s %s' \
	       % ( binPath, object_param['isoResolution'], object_param['isoNlevel'], object_param['isoMargin'], \
	       	   object_param['isoAlpha'], object_param['isoBeta'], object_param['objFile'], object_param['tetFile'] );
	print cmd;
	os.system( '%s > %s' % ( cmd, redirect ) );
	print("=== Done ===")

	# Step two - generate the mass and stiffness matrix
	print("=== Generating the mass and stiffness matrix ===")
	cmd = '%s/elasticity_solver %s %f %f' \
	       % ( binPath, object_param['tetFile'], material['youngModulus'], material['poissonRatio'] );
	print cmd;
	os.system( '%s > %s' % ( cmd, redirect ) );
	print("=== Done ===")

	# This step also generates plate.tet.geo.txt, but we want to rename this
	os.system( 'mv -f %s.geo.txt %s.geo.txt' % ( object_param['tetFile'], object_param['objName'] ) );

	# Step three - linear modal analysis
	#
	# Using a threshold of 1.0 here seems to work
	print("=== Performing linear modal analysis ===")
	cmd = ("%s/arpack-eigensolver -n %d -t 1.0 -s %s_stiffness.mat -m %s_mass.mat "
	       "-o %s.modes -v") % ( binPath, object_param['numEigs'], object_param['tetFile'], object_param['tetFile'],\
	       						 object_param['objName'] );
	print cmd;
	os.system( '%s' % ( cmd) );
	print("=== Done ===")

	# Step four - generate input files for the Helmholtz radiation problems
	os.system( 'mkdir -p fbem_in' );
	cmd = '%s/fbem_input_gen -m %s.modes -d %f -t %s -o fbem_in/input-%s_please_d.txt' \
	       % ( binPath, object_param['objName'], material['density'], object_param['tetFile'], r"%d");
	print cmd;
	os.system( '%s > %s' % ( cmd, redirect ) );


	# Step five - write the object center of mass, rigid inertia matrix, and
	# material parameters to files that will be needed in other parts of the
	# pipeline
	print("=== Writing rigid data ===")
	cmd = '%s/write-mass-center %s %s' % ( binPath, object_param['tetFile'], object_param['dataPrefix'] );
	print cmd;
	os.system( '%s > %s' % ( cmd, redirect ) );

	cmd = '%s/write-inertia %s %s' % ( binPath, object_param['tetFile'], object_param['dataPrefix'] );
	print cmd;
	os.system( '%s > %s' % ( cmd, redirect ) );

	# This part assumes that octave is installed.  This command could be
	# changed somewhat to work with Matlab
	cmd = ("octave --eval \"addpath('%s'); "
	       "make_material_parameters(%f, %f, '%s_materialProperties.vector');\"") \
	       % ( matlabPath, material['youngModulus'], material['poissonRatio'], object_param['dataPrefix'] );
	print cmd;
	os.system( '%s > %s' % ( cmd, redirect ) );
	print("=== Done ===")

	# # Step six - wrap everything up in a tarball so that it can be conveniently
	# # moved to the hydra cluster.
	# #
	# # As a final step, move this file to /hydra/S1/impact_snd/fastbem/models and
	# # untar it there (e.g., tar xzvf plate.tar.gz).  This should create a 'plate'
	# # directory.  Inside this directory, you can run the following command:
	# #     ./../../scripts/fbem_solve.pl . <modeID>
	# # This will solve the boundary element problem for the given mode ID, placing
	# # the results in fbem_ret
	# # 
	# cmd = ("tar --transform 's,^,%s/,' -czvf %s.tar.gz fbem_in/ "
	#        "freqs.txt %s.geo.txt") % ( objName, objName, objName );
	# print cmd;
	# os.system( '%s > %s' % ( cmd, redirect ) );
