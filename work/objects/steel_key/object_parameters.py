import os

def object_parameters():
	return {
		'steel_key':{	
				'objName': 'steel_key',
				'objPath': os.path.dirname(__file__),
				'objFile': os.path.abspath('../../meshes/key/key.obj'),
				'tetFile': os.path.dirname(__file__)+'/steel_key.tet',
				'dataPrefix': os.path.dirname(__file__)+'/steel_keyRigid',
				'material': 'steel',
				'numEigs': 50,
				'isoResolution': 7,
				'isoNlevel': 3,
				'isoMargin': 7,
				'isoAlpha': 0.25,
				'isoBeta': 0.42978,
				'wave_params':	{
									"field_resolution": 100,
									"grid_resolution": 100,
									"grid_scale": 2.0,
									"time_step_freq": 440000,
									"substeps": 8,
									"radius_multipole": 1.5
								}
		}
	}