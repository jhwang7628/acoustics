#!/usr/bin/env python
import sys
import os
################################################################################
## This code create tet mesh from obj
################################################################################

################################################################################
################################################################################
def ReadMaterialFile(filename): 
    lines = open(filename, 'r').readlines()
    for l in lines: 
        tokens = l.split()
        if tokens[0] == '#':
            if tokens[1] == 'MATERIAL_PROPERTIES': 
                name = tokens[-1]
        else: 
            rho = float(tokens[0])
            youngs = float(tokens[1])
            poisson = float(tokens[2])
    print 'Parsed material for file: %s' %(tokens[-1])
    print ' density        : %f' %(rho)
    print ' youngs modulus : %e' %(youngs)
    print ' poisson ratio  : %.2f' %(poisson)
    return rho, youngs, poisson

################################################################################
################################################################################
if __name__ == '__main__': 
    if len(sys.argv) != 4: 
        print '**Usage: %s <obj_prefix> <material_file> <num_eigenvalues>' %(sys.argv[0])
        sys.exit()
    
    ## User defined settings
    installPath = '/home/jui-hsien/code/acoustics'
    binPath = '%s/build_release/bin' %(installPath)
    redirect = '/dev/stdout';
               
    ## Automatic 
    print '################################################################################'
    print '## Initialization'
    print '################################################################################'
    objName = sys.argv[1]
    matName = sys.argv[2]
    numEigs = int(sys.argv[3])
    tetFile = '%s.tet' %(objName)
    bin_elasticity='%s/elasticity_solver' %(binPath)
    bin_arpack_eigensolver='%s/arpack-eigensolver' %(binPath)
    density, youngsModulus, poissonRatio = ReadMaterialFile(matName)
    if not os.path.isfile(bin_elasticity) or not os.path.isfile(bin_arpack_eigensolver): 
        print '**ERROR** some binaries required are not found in path %s' %(binPath)
    print '\n\n'
    
    ## Pipeline starts
    print '################################################################################'
    print '## Generate mass and stiff matrix'
    print '################################################################################'
    cmd = '%s %s %f %f' %(bin_elasticity, tetFile, youngsModulus, poissonRatio);
    print cmd;
    print '........................................'
    os.system('%s > %s' %(cmd, redirect));
    print '\n\n'
    
    # This step also generates <objName>.tet.geo.txt, but we want to rename this
    os.system('mv -f %s.geo.txt %s.geo.txt' %(tetFile, objName));
    
    print '################################################################################'
    print '## Linear modal analysis'
    print '################################################################################'
    # Using a threshold of 1.0 here seems to work
    cmd = ("%s -n %d -t 1.0 -s %s_stiffness.mat -m %s_mass.mat -o %s.modes -v") %(bin_arpack_eigensolver, numEigs, tetFile, tetFile, objName);
    print cmd;
    print '........................................'
    os.system('%s > %s' %(cmd, redirect));
    print '\n\n'
