import set_lib_path

import pipeline

import scene_parameters

if __name__ == '__main__':
	objects = scene_parameters.objects()
	objDict = scene_parameters.objectDict()
	sim_parameters = scene_parameters.sim_parameters()
	pipeline.run_simulation(sim_parameters, objDict, objects)

		