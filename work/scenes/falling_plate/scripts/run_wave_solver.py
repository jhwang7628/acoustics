import set_lib_path
from pylab import *
import pipeline
import math
import scene_parameters

def absVec(vec):
	return math.sqrt(sum([i*i for i in vec]))

def pre_wave_strategy(solv, sim_params, obj, mode, frequency):
	print "HI!"
	print "Radius %f"%(solv.radius)

def post_wave_strategy(solv, sim_params, obj, mode, frequency):
	print frequency
	nbar = max(5, 0.22*2*math.pi*frequency/(2.0*343.0))
	nbar = min(nbar, 20)
	nbar = int(math.floor(nbar))
	print "nbar = %d" % (nbar)
	solv.nbar = nbar

	solv.setEndTime(600)
	print "Radius %f"%(solv.radius)

import matplotlib.pyplot as plt
def plotdata(data, filename):
	maxi = max([max([abs(d) for d in item]) for item in data])
	plt.clf()
	plt.imshow(data, aspect='auto', cmap=get_cmap('seismic'), origin='lower', vmin=-maxi, vmax=maxi)
	plt.colorbar()
	plt.savefig(filename)

def plot_wave(solv, sim_params, obj, moded, frequency):
	size = solv.cellDivisions()
	z = size[2]/2

	data = []
	for x in range(0, size[0]):
		line = []
		for y in range(0, size[1]):
			cellData = solv.cellData(x, y, z)
			gradient = solv.gradientAt(x, y, z)
			datum = absVec(gradient)*frequency*2*math.pi
			if not cellData[2]:
				datum = 0
			line.append(datum)
		data.append(line)
	plotdata(data, "plot_%d_z_slice_grad.png"%(moded))

	y = size[1]/2

	data = []
	for x in range(0, size[0]):
		line = []
		for z in range(0, size[2]):
			cellData = solv.cellData(x, y, z)
			gradient = solv.gradientAt(x, y, z)
			datum = absVec(gradient)*frequency*2*math.pi
			if not cellData[2]:
				datum = 0
			line.append(datum)
		data.append(line)
	plotdata(data, "plot_%d_y_slice_grad.png"%(moded))

	x = size[0]/2

	data = []
	for y in range(0, size[1]):
		line = []
		for z in range(0, size[2]):
			cellData = solv.cellData(x, y, z)
			gradient = solv.gradientAt(x, y, z)
			datum = absVec(gradient)*frequency*2*math.pi
			if not cellData[2]:
				datum = 0
			line.append(datum)
		data.append(line)
	plotdata(data, "plot_%d_x_slice_grad.png"%(moded))

	y = size[1]/2
	data = []
	for x in range(0, size[0]):
		line = []
		for z in range(0, size[2]):
			cellData = solv.cellData(x, y, z)
			datum = cellData[0]
			line.append(datum)
		data.append(line)
	plotdata(data, "plot_%d_y_slice_temp.png"%(moded))

	x = size[0]/2
	data = []
	for y in range(0, size[1]):
		line = []
		for z in range(0, size[2]):
			cellData = solv.cellData(x, y, z)
			datum = cellData[0]
			line.append(datum)
		data.append(line)
	plotdata(data, "plot_%d_x_slice_temp.png"%(moded))

	z = size[2]/2
	data = []
	for x in range(0, size[0]):
		line = []
		for y in range(0, size[2]):
			cellData = solv.cellData(x, y, z)
			datum = cellData[0]
			line.append(datum)
		data.append(line)
	plotdata(data, "plot_%d_z_slice_temp.png"%(moded))



if __name__ == '__main__':
	objDict = scene_parameters.objectDict()
	for objType in objDict:
		print objType
		pipeline.compute_multipole_coefficients(scene_parameters.sim_parameters(),\
												objDict[objType],\
												pre_wave_strategy, post_wave_strategy, plot_wave)

		
