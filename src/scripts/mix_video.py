#!/usr/bin/python
#
# Estimates impact timescales using Hertz contact theory

import sys
import os

def mix_video(sim_params):
	cmd = 'mencoder no_audio_%s -o %s -ovc copy -oac copy -audiofile %s' % (sim_params['videoFile'], sim_params['videoFile'], sim_params['audioFile'])

	print( "=== Mixing" );
	os.system(cmd);
	print( "=== Done" );