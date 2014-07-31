import multipole
import solver
import os
import math
import mode
import time
import acoustic_templates
from material_parameters import materials

def null_wave_strategy(sim_params, obj, mode, frequency):
	print "<WARNING>: Null Wave Strategy is being used!"

def compute_multipole_coefficients(sim_params, obj, wave_strategy_pre_sdf=null_wave_strategy, wave_strategy_pre_init=null_wave_strategy):
	timing = {}
	timing["average"] = 0
	timing["min"] = 100000
	timing["max"] = 0
	timing["total"] = 0


	configfile = obj['objName']+"__wavegen_config.xml"
	with open(configfile, 'w+') as conffile:
		conffile.write(acoustic_templates.wave_config(sim_params, obj))


	modefile = obj['objName']+".modes"
	material = materials()[objct['material']]

	modeData = mode.ModeData()
	modeData.read(modefile)

	nmodes = modeData.numModes()

	files = multipole.StringVector()

	for mod in range(0, nmodes):
		tic = time.time()
		print ">>>>>> DOING mode " + str(mod)
		frequency = math.sqrt(modeData.omegaSquared[mod]/material['density'])/(2*math.pi)
		print "Mode " + str(mod) + ": " + str(frequency) + " Hz"

		if frequency < sim_params['minFreq'] or frequency > sim_params['maxFreq']:
			print "><> Ignore mode %d" % (mod)
			continue

		solv = solver.PAT_Solver(configfile)
		solv.mode = mod
		# Someone has to define a strategy
		wave_strategy_pre_sdf(sim_params, obj, mod, frequency)
		solv.initSDF()
		# Someone has to define a strategy
		wave_strategy_pre_init(sim_params, obj, mode, frequency)
		solv.initSolver()
		solv.runSolver()
		fil = "%s-mode-%d.multipole"%(objName[objName], mod)
		solv.saveToFile(fil)
		files.append(fil)

		del solv

		print toc
		timing["total"] += toc
		timing["min"] = min(timing["min"], toc)
		timing["max"] = max(timing["max"], toc)

	multidata = multipole.mergeMultipoleModeDataFiles(files, 343, 0, 0, 0)
	multipole.saveMultipoleDataToFile("%s.multipoleAll"%(obj['objName']), multidata)

	# Cleaning data
	for f in files:
		os.remove(f)

	print timing