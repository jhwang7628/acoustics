import os
import matplotlib.pyplot as plt
import sys

def main(dir):

	os.chdir(dir)
	data = {}

	for tfile in os.listdir('.'):
		if 'tfv' in tfile and 'alltransfer' not in tfile:
			tdata = []
			with open(tfile) as transfer:
				for line in transfer:
					w = line.split()
					tdata.append(complex(float(w[0]), float(w[1])))

			data[tfile] = tdata
	for tdata in data:
		toplot = [abs(transfer) for transfer in data[tdata]]
		plt.semilogy(toplot)
		plt.ylabel(dir)

	plt.show()
	
if __name__ == '__main__':
	main(sys.argv[1])