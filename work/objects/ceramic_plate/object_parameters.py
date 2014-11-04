import os

def object_parameters():
	return {
		'ceramic_plate':{	
				'objName': 'ceramic_plate',
				'objPath': os.path.dirname(__file__),
				'objFile': os.path.abspath('../../meshes/plate/plate.obj'),
				'tetFile': os.path.dirname(__file__)+'/ceramic_plate.tet',
				'dataPrefix': os.path.dirname(__file__)+'/ceramic_plateRigid',
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
									"grid_scale": 3.2,
									"time_step_freq": 192000,
									"substeps": 8,
									"radius_multipole": 1.8
								}
		}
	}