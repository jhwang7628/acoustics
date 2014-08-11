import os, sys

solver = '/home/rf356/FastBEM/program/FastBEM_Acoustics.exe'

def main():
	tests = [o for o in os.listdir('.') if os.path.isdir(o)]
	now = os.getcwd()
	print now
	for test in tests:
		os.chdir(now)
		dapath = test
		modes = [m for m in os.listdir(dapath)]
		for mode in modes:
			os.chdir(now)
			damode = os.path.join(dapath, mode)
			damode = os.path.abspath(damode)
			os.chdir(damode)
			cmd = "echo . | %s " % (solver)
			print cmd
			os.system(cmd)
			

if __name__ == '__main__':
	main()