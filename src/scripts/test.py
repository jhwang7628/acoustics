import multipole
import solver
import os
import math
import mode
import time
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


os.chdir("../src/work/plate")

data = {}
data["dir"] = "."
data["configfile"] = "default.xml"
data["modefile"] = "plate.modes"

print data

os.listdir(data["dir"])

modeData = mode.ModeData()
modeData.read(data["modefile"])

modes = 0

tempsolv = solver.PAT_Solver(data["configfile"])
density = tempsolv.density
del tempsolv

for mod in range(0, modeData.numModes()):
	frequency = math.sqrt(modeData.omegaSquared[mod]/density)/(2*math.pi)
	if frequency > 20000:
		break
	modes = modes + 1

for mod in range(0, modes):
	tic = time.time()

	print ">>>>>> DOING mod " + str(mod)
	solv = solver.PAT_Solver(data["configfile"])
	frequency = math.sqrt(modeData.omegaSquared[mod]/solv.density)/(2*math.pi)

	print "Mode " + str(mod) + ": " + str(frequency) + " Hz"
	solv.initSDF()
	solv.mode = mod
	k = 2*math.pi*frequency/solv.wave_speed
	sugnbar = max(4, int(k*solv.radius/2))
	sugnbar = 5
	solv.nbar = sugnbar
	solv.setEndTime(500)
	solv.initSolver()
	if frequency < 14000:
		solv.runSolver()
	solv.saveToFile("pplate"+str(mod)+".multipole")
	del solv
	print "<<<<<< FINISHED mod" + str(mod)

	toc = time.time() - tic

	print toc
	timing["total"] += toc
	timing["min"] = min(timing["min"], toc)
	timing["max"] = max(timing["max"], toc)

print "USED " + str(modes) + " modes"

timing["average"] = timing["total"]/modes

files = multipole.StringVector()
f = [f for f in os.listdir(data["dir"])
		if f.startswith("pplate")]
files.extend(f)
print f
multidata = multipole.mergeMultipoleModeDataFiles(files, 343, 0, 0, 0)
multipole.saveMultipoleDataToFile("plate.multipoleAll", multidata)

print timing