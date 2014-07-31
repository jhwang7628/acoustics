import sys
import multipole

def evaluate_transfer(sim_params, objectDict, objects):
	os.system( 'mkdir -p fbem_ret' );
	for i in range(0, len(objects)):
		obj = objects[i]
		inputFile = sim_params['outputFilePattern']%(i)
		inputMultipoleFile = obj['object']+".multipoleAll"
		multipoleData = multipole.loadMultipoleDataFromFile(inputMultipoleFile)
		numModes = len(multipoleData.modes)
		tModes = max(numModes, minimumModes)
		for mode in range(0, tmodes):
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
							val = multipoleData.estimateModeAt(mode, x, y, z)
					
					outData.append(val)

			fname = sim_params['transferFilePattern']%(i, mode)

			with open(fname, 'w+') as outputData:
				for data in outData:
					outputData.write(str(data.real) + " " + str(data.imag) + "\n")