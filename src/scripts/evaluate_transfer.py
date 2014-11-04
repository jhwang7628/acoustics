import os, sys
import multipole

def evaluate_transfer(sim_params, objectDict, objects):
	os.system( 'mkdir -p fbem_ret' );
	for i in range(0, len(objects)):
		obj = objects[i]
		inputFile = sim_params['outputFilePattern']%(i)
		inputMultipoleFile = obj['object']+".multipoleAll"
		multipoleData = multipole.loadMultipoleDataFromFile(inputMultipoleFile)
		numModes = len(multipoleData.modes)
		tModes = max(numModes, sim_params['minModes'])
		center = (multipoleData.cx, multipoleData.cy, multipoleData.cz)
		print center
		for mode in range(0, tModes):
			outData = []
			if mode in multipoleData.modes:
				frequency = multipoleData.modes[mode].frequency()
				if frequency > 20000:
					break

			with open(inputFile) as fieldPoints:
				for line in fieldPoints:
					w = line.split()
					x = float(w[1])
					y = float(w[2])
					z = float(w[3])

					val = 0 + 0j
					if mode in multipoleData.modes:
						frequency = multipoleData.modes[mode].frequency()
						if frequency >= sim_params['minFreq'] and frequency <= sim_params['maxFreq']:
							val = multipoleData.estimateModeAt(mode, x+center[0], y+center[1], z+center[2])
					
					# if frequency > 8000:
					# 	val = val*1000

					outData.append(val)

			fname = sim_params['transferFilePattern']%(i, mode)

			with open(fname, 'w+') as outputData:
				for data in outData:
					outputData.write(str(data.real) + " " + str(data.imag) + "\n")