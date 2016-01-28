import os

def object_parameters():
	return {
		'ceramic_mug':{	
				'objName': 'ceramic_mug',
				'objPath': os.path.dirname(__file__),
				'objFile': os.path.abspath('../../meshes/mug/mug.obj'),
				'tetFile': os.path.dirname(__file__)+'/ceramic_mug.tet',
				'dataPrefix': os.path.dirname(__file__)+'/ceramic_mugRigid',
				'material': 'ceramic',
				'numEigs': 50,
				'isoResolution': 6,
				'isoNlevel': 3,
				'isoMargin': 7,
				'isoAlpha': 0.25,
				'isoBeta': 0.42978,
				'wave_params':	{
									"field_resolution": 150,
									"grid_resolution": 200,
									"grid_scale": 2.0,
									"time_step_freq": 440000,
									"substeps": 8,
									"radius_multipole": 1.5
								}
		}
	}