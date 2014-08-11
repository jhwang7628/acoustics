import set_lib_path
import pipeline
import object_parameters

if __name__ == '__main__':
	obj_d = object_parameters.object_parameters()
	for obj in obj_d:
		print obj
		pipeline.precompute_modes(obj_d[obj])