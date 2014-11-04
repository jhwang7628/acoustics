import os

def templ(k):
	t = """Parameters for fast multipole method:
100 	10      %d	1.E-4	800000000    3500 ! maxl     levmx    nexp  tolerance  nstk   iteration
60000	50000	100	1	%d                ! ncellmx  nleafmx  maxk  ifmm       nexpf
60000	50000	100	0	10                ! ncellmx  nleafmx  maxk  ifmm       nexpf
>>>>>>> tim-gpu


-----------------------------------------------------------------------------------------

Definitions of the above parameters:

maxl: 	    maximum number of elements in a leaf
levmx:	    maximum number of tree levels
nexp:	    order of fast multipole expansions (p)
tolerance:  tolerance for convergence used in the iterative solver
nstk:       size of the array used to store coefficients in preconditioner (value > 0 or = 0)
iteration:  maximum number of iterations allowed for the iterative solver
ncellmx:    maximum number of cells allowed in the tree
nleafmx:    maximum number of leaves allowed in the tree
maxk:       maximum dimension of Krylov subspace used in the iterative solver
ifmm:	    apply FMM in field evaluation (0-no; 1-yes)
nexpf:	    number of expansion terms in field evaluation if FMM is used


-----------------------------------------------------------------------------------------

Default Values of these parameters:

100	10	6	1.E-4	100000000    500  ! maxl     levmx    nexp  tolerance  nstk   iteration
60000	50000	100	0	10                ! ncellmx  nleafmx  maxk  ifmm       nexpf


-----------------------------------------------------------------------------------------
"""
	return t%(k, k)

import math

def main():
	pointsFileTemp = 'scenes/%s/field_pts/fp-0.txt'
	inputFileTemp = 'objects/%s/fbem_in/input-%d.txt'
	outputDirTemp = 'fbem_test/%s/mode_%d'

	dalist = [
				("ceramic_mug", "falling_mug", 28, 0.20),
				("ceramic_plate", "falling_plate", 38, 0.22),
				("steel_key", "falling_key", 3, 0.08)	
			 ]

	for item in dalist:
		print item
		points = []

		pointsFile = pointsFileTemp % (item[1])

		with open(pointsFile) as ptsFile:
			for line in ptsFile:
				w = line.split()
				points.append((w[1], w[2], w[3]))

		for mode in range(0, item[2]+1):
			outputLines = []
			inputFile = inputFileTemp % (item[0], mode)

			frequency = 1000.0

			with open(inputFile) as iptFile:
				for (i, line) in enumerate(iptFile):
					if i == 3:
						w = line.split()
						outputLines.append("%s %s %d %s\n" % (w[0], w[1], len(points), w[2]))
					else:
						outputLines.append(line)
					
					if i == 7:
						w = line.split()
						frequency = float(w[1])

					if "Field Points" in line:
						for i, point in enumerate(points):
							outputLines.append("%d %s %s %s\n" % (i+1, point[0], point[1], point[2]))

			print "%s - Mode %d of %d - %f Hz" % (item[0], mode, item[2], frequency)
			outputDir = outputDirTemp % (item[0], mode)

			if not os.path.exists(outputDir):
				os.makedirs(outputDir)

			with open(outputDir + "/input.dat", 'w+') as out:
				for line in outputLines:
					out.write(line)

			k = max(5, item[3]*2*math.pi*frequency/(4.0*343.0))
			k = int(math.floor(k))
			with open(outputDir + "/input.fmm", 'w+') as out:
				out.write(templ(k))

if __name__ == '__main__':
	main()
