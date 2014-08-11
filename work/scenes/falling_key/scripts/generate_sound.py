import set_lib_path

import pipeline

import scene_parameters

if __name__ == '__main__':
	objects = scene_parameters.objects()
	objDict = scene_parameters.objectDict()
	sim_parameters = scene_parameters.sim_parameters()
	pipeline.prep_transfer_eval(sim_parameters, objDict, objects)
	pipeline.evaluate_transfer(sim_parameters, objDict, objects)
	pipeline.estimate_timescales(sim_parameters, objDict, objects)
	pipeline.generate_sound(sim_parameters, objDict, objects)
	pipeline.mix_video(sim_parameters)


		