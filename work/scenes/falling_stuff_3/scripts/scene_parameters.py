import set_lib_path
from objects import ceramic_mug
from objects import steel_key
from objects import ceramic_plate

def objectDict():
	objDict = {}
	objDict.update(ceramic_mug.object_parameters())
	objDict.update(ceramic_plate.object_parameters())
	objDict.update(steel_key.object_parameters())
	return objDict

def objects():
	return [
		{
			'name': 'mug_0',
			'object': 'ceramic_mug',
			'position': (0.0, 0.402, 0.0),
			'velocity': (0.0, 0, 0),
			'rotation': (1, 0, 0, 0)
		},
		{
			'name': 'mug_1',
			'object': 'ceramic_mug',
			'position': (0.0, 0.802, 0.0),
			'velocity': (0.0, 0.0, 0),
			'rotation': (1, 0, 0, 0)
		},
		{
			'name': 'key_0',
			'object': 'steel_key',
			'position': (0.0, 0.602, 0.0),
			'velocity': (0.0, 0, 0),
			'rotation': (1, 0, 0, 0)
		},
		{
			'name': 'key_1',
			'object': 'steel_key',
			'position': (0.0, 1.002, 0.0),
			'velocity': (0.0, 0.0, 0),
			'rotation': (1, 0, 0, 0)
		},
		{
			'name': 'plate_1',
			'object': 'ceramic_plate',
			'position': (-0.1, 0.2, 0.0),
			'velocity': (0.0, 0.0, 0),
			'rotation': (-0.7, 0.7, 0, 0)
		}
	]

def sim_parameters():
	return { 
				'simulationFile': 'default.cfg',
				'soundConfigFile': 'example.cfg',
				'displaceFile': 'displace.bin',
				'frameDir': 'frames',
				'videoFile': 'mug.mpg',
				'audioFile': 'mug.wav',
				'massCenterFile': 'mc.txt',
				'mapFile': 'obj_model_map.txt',
				'outputFilePattern': 'field_pts/fp-%d.txt',
				'transferFilePattern': 'fbem_ret/tfv-%d_%d.txt',
				'sound_rate': 192000,
				'sim_step': 0.001,
				'simulation_length': 2.5,
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
