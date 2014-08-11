import os, sys

installPath = '/home/rf356/Desktop/Code/acoustics';

libpath = os.path.abspath('%s/gcc-build/lib'%(installPath))
workpath = os.path.abspath('%s/work/'%(installPath))

sys.path.append(libpath)
sys.path.append(workpath)