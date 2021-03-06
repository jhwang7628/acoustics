import set_lib_path

import pipeline

import scene_parameters

def pre_wave_strategy(solv, sim_params, obj, mode, frequency):
	print "HI!"

def post_wave_strategy(solv, sim_params, obj, mode, frequency):
	print "nbar = 5"
	solv.nbar = 5
	solv.setEndTime(500)

if __name__ == '__main__':
	objDict = scene_parameters.objectDict()
	for objType in objDict:
		print objType
		pipeline.compute_multipole_coefficients(scene_parameters.sim_parameters(),\
												objDict[objType],\
												pre_wave_strategy, post_wave_strategy)

		