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
	print "Radius %f"%(solv.radius)
	solv.setEndTime(800)



import matplotlib.pyplot as plt
from pylab import *

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
			datum = gradient[2]
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
			datum = gradient[2]
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
			datum = gradient[2]
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

		
