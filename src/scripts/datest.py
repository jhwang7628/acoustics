import os, sys

installPath = '/home/rf356/Desktop/Code/acoustics';
libpath = os.path.abspath('%s/gcc-build/lib'%(installPath))
sys.path.append(libpath)
import multipole
from material_parameters import materials
from params import parameters
print materials()
print parameters()