import re
import os

def extract_number(s,notfound='NOT_FOUND'):
    regex=r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
    return re.findall(regex,s)

def parse_transfer(inputFile, outputFile):
	data = []
	dafile = open(inputFile)
	try:
		state = 0
		for line in dafile:
			if state == 1:
				w = extract_number(line)
				if len(w) == 7:
					data.append((float(w[1]), float(w[2])))

			if "Point #" in line:
				state = 1
			if "Sound power" in line:
				state = 0
	finally:
		dafile.close()

	dafile = open(outputFile, 'w+')
	try:
		for datum in data:
			dafile.write(str(datum[0]) + " " + str(datum[1]) + "\n")
	finally:
		dafile.close()

def main():

	tests = [o for o in os.listdir('.') if os.path.isdir(o)]
	now = os.getcwd()
	print now
	
	if not os.path.exists("transfer"):
		os.makedirs("transfer")

	for test in tests:
		if "transfer" in test:
			continue
			
		os.chdir(now)
		dapath = test
		modes = [m for m in os.listdir(dapath)]
		outPath = os.path.join("transfer", test)
		if not os.path.exists(outPath):
			os.makedirs(outPath)

		for mode in modes:
			os.chdir(now)
			modeNumber = mode[5:]
			damode = os.path.join(dapath, mode)
			dafile = os.path.join(damode, "output_result.dat")

			outputFile = "tfv-0_%s.txt" % (modeNumber)
			daOutFile = os.path.join(outPath, outputFile)
			print "Test %s - Mode %s" % (test, modeNumber)
			parse_transfer(dafile, daOutFile)

if __name__ == '__main__':
	main()


