def sim_parameters():
	return { 
				'simulationFile': 'default.cfg',
				'soundConfigFile': 'example.cfg',
				'displaceFile': 'displace.bin',
				'frameDir': 'frames',
				'videoFile': 'plate.mpg',
				'audioFile': 'plate.wav',
				'massCenterFile': 'mc.txt',
				'mapFile': 'obj_model_map.txt',
				'outputFilePattern': 'field_pts/fp-%d.txt',
				'transferFilePattern': 'fbem_ret/tfv-%d_%d.txt',
				'sound_rate': 192000,
				'sim_step': 0.001,
				'simulation_length': 1.1,
				'frameFrequency': 10.0,
				'transferFrequency': 1000.0,
				'listeningPosition': (1.0, 1.0, 1.0),
				'forceScalePrefix': 'forceScales',
				'impulseFile': 'impulses.txt',
				'modalImpulse': 'modalImpulses.txt',
				'minFreq': 300.0,
				'maxFreq': 14000.0,
				'minModes': 59
			}


def objectDict():
	return {
		'ceramic_plate':{	
			'objName': 'ceramic_plate',
			'objFile': 'plate.obj',
			'tetFile': 'ceramic_plate.tet',
			'dataPrefix': 'ceramic_plateRigid',
			'material': 'ceramic',
			'numEigs': 100,
			'isoResolution': 7,
			'isoNlevel': 3,
			'isoMargin': 7,
			'isoAlpha': 0.25,
			'isoBeta': 0.42978,
			'wave_params':	{
								"field_resolution": 150,
								"grid_resolution": 300,
								"grid_scale": 2.95,
								"time_step_freq": 192000,
								"substeps": 8,
								"radius_multipole": 4.0
							}
		}
	}

def objects():
	return [
		{
			'name': 'plate_0',
			'object': 'ceramic_plate',
			'position': (0.0, 0.602, 0.0),
			'rz' : 90.0,
			'rotation': (0.02895, -0.267, 0.337, -0.902)
		}
	]

