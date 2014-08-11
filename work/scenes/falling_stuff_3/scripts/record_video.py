import set_lib_path

import pipeline

import scene_parameters

if __name__ == '__main__':
	objects = scene_parameters.objects()
	objDict = scene_parameters.objectDict()
	sim_parameters = scene_parameters.sim_parameters()
	pipeline.make_frames(sim_parameters)
	pipeline.encode_frames(sim_parameters)

		