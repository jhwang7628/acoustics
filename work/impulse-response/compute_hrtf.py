import argparse
import os
import sys
from subprocess import call
# 0 -> Left
# 1 -> Right

def convert_to_earfile(listener, ear):

        listener_=os.path.split(listener)

	fin = open(listener, 'r')
	lines = fin.readlines()	
	fin.close()

	nfname = ""

        sourceL = "%lf %lf %lf\n" %(0.0, 0.0,  0.08)
        sourceR = "%lf %lf %lf\n" %(0.0, 0.0, -0.08)

        # sourceL = "%lf %lf %lf\n" %(0.0, 0.0, 0.0)
        # sourceR = "%lf %lf %lf\n" %(0.0, 0.0, 1.25)


        sout = open(listener_[0]+'/source_position.txt', 'w')
        sout.write('----\n')

	if ear == 0:
		nfname = "temp_left_"
		lines[0] = lines[0]%("left_%s")
		lines.insert(2, sourceL)
                sout.write(sourceL)
	else:
		nfname = "temp_right_"
		lines[0] = lines[0]%("right_%s")
		lines.insert(2, sourceR)
                sout.write(sourceR) 

	nfname = listener_[0] + '/' + nfname + listener_[1]
	fout = open(nfname, 'w')
	fout.writelines(lines)
	fout.close()

        sout.close()

	return nfname

def main(args):

        outdir = open(args.listener_file,'r').readline().split('/')[0]

        call( 'cp %s %s %s' %(args.wave_file, args.listener_file, outdir), shell=True )

	if args.visual_ear:
                os.environ["PATH"] += os.pathsep + ("/home/jui-hsien/code/acoustics/build_ui/bin")
		if args.visual_ear == "left":
		    efile = convert_to_earfile(args.listener_file, 0)
		    os.system("precompute-impulse-response %s %s"%(args.wave_file, efile))
		elif args.visual_ear == "right":
		    efile = convert_to_earfile(args.listener_file, 1)
		    os.system("precompute-impulse-response %s %s"%(args.wave_file, efile))
	else:
                os.environ["PATH"] += os.pathsep + ("/home/jui-hsien/code/acoustics/build/bin")
		# print "Left Ear"
		# efile = convert_to_earfile(args.listener_file, 0)
		# os.system("precompute-impulse-response-text %s %s"%(args.wave_file, efile))

		print "Right Ear"
		efile = convert_to_earfile(args.listener_file, 1)
		os.system("precompute-impulse-response-text %s %s"%(args.wave_file, efile))


        call( 'rm temp*', shell=True )

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Run GPU wave solver')
	parser.add_argument("wave_file", help='setup wavesolver')
	parser.add_argument("listener_file", help='setup master listening positions')
        parser.add_argument("--centroid_file", help='cell centroid files')
        parser.add_argument("--centroid_file2", help='cell centroid files 2')
	parser.add_argument("--visual_ear", help='use GUI interface to visualize')
	args = parser.parse_args()

	main(args)
