import argparse
import os
# 0 -> Left
# 1 -> Right

def convert_to_earfile(listener, ear):
	fin = open(listener, 'r')
	lines = fin.readlines()	
	fin.close()

	nfname = ""

	if ear == 0:
		nfname = "temp_left_"
		lines[0] = lines[0]%("left_%s")
		lines.insert(2, "%lf %lf %lf\n"%(-0.0763, 0.003, -0.005))
	else:
		nfname = "temp_right_"
		lines[0] = lines[0]%("right_%s")
		lines.insert(2, "%lf %lf %lf\n"%(0.069, 0.0, 0.0))

	nfname = nfname + listener
	fout = open(nfname, 'w')
	fout.writelines(lines)
	fout.close()

	return nfname

def main(args):

	if args.visual_ear:
		if args.visual_ear == "left":
			efile = convert_to_earfile(args.listener_file, 0)
			os.system("../../gcc-build2/bin/precompute-impulse-response %s %s"%(args.wave_file, efile))
		elif args.visual_ear == "right":
			efile = convert_to_earfile(args.listener_file, 1)
			os.system("../../gcc-build2/bin/precompute-impulse-response %s %s"%(args.wave_file, efile))
	else:
		print "Left Ear"
		efile = convert_to_earfile(args.listener_file, 0)
		os.system("../../gcc-build2/bin/precompute-impulse-response-text %s %s"%(args.wave_file, efile))

		# print "Right Ear"
		# efile = convert_to_earfile(args.listener_file, 0)
		# os.system("../../gcc-build/bin/precompute-impulse-response-text %s %s"%(args.wave_file, efile))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("wave_file")
	parser.add_argument("listener_file")
	parser.add_argument("--visual_ear")
	args = parser.parse_args()

	main(args)
