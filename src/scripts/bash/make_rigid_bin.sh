#!/bin/bash 
# Just saving myself some troubles on building individual targets for modal stuff

TARGET="rigidsim rigidsim-gui write-mass-center"

# make -j40 isostuffer elasticity_solver arpack-eigensolver
for p in ${TARGET}; do 
    echo -e "\nMaking target $p" 
    make -j40 $p
done 
