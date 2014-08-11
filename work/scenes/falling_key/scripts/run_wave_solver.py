import set_lib_path

import pipeline

import scene_parameters
import math

def pre_wave_strategy(solv, sim_params, obj, mode, frequency):
	print "HI!"
	print "Radius %f"%(solv.radius)

def post_wave_strategy(solv, sim_params, obj, mode, frequency):
	if frequency < 7000:
		solv.nbar = 5
	else:
		solv.nbar = max(5, int(math.floor(2*(solv.radius*solv.scaleRadius)*(frequency*math.pi)/343.0)))
	print "nbar = %d" % (solv.nbar)
	solv.setEndTime(800)
	print "Radius %f"%(solv.radius)

if __name__ == '__main__':
	objDict = scene_parameters.objectDict()
	for objType in objDict:
		print objType
		pipeline.compute_multipole_coefficients(scene_parameters.sim_parameters(),\
												objDict[objType],\
												pre_wave_strategy, post_wave_strategy)

		
