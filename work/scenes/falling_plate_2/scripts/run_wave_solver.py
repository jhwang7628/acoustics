import set_lib_path

import pipeline
import math
import scene_parameters

def pre_wave_strategy(solv, sim_params, obj, mode, frequency):
	print "HI!"
	print "Radius %f"%(solv.radius)

def post_wave_strategy(solv, sim_params, obj, mode, frequency):
	print frequency
	nbar = max(5, 0.22*2*math.pi*frequency/(4.0*343.0))
	nbar = int(math.floor(nbar))
	print "nbar = %d" % (nbar)
	solv.nbar = nbar

	solv.setEndTime(600)
	print "Radius %f"%(solv.radius)

if __name__ == '__main__':
	objDict = scene_parameters.objectDict()
	for objType in objDict:
		print objType
		pipeline.compute_multipole_coefficients(scene_parameters.sim_parameters(),\
												objDict[objType],\
												pre_wave_strategy, post_wave_strategy)

		
