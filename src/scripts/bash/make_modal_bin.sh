#!/bin/bash 
# Just saving myself some troubles on building individual targets for modal stuff

TARGET="isostuffer elasticity_solver arpack-eigensolver"
if [ $# -eq 1 ]; then 
    TARGET="${TARGET} fbem_input_gen write-mass-center write-inertia"
fi 

# make -j40 isostuffer elasticity_solver arpack-eigensolver
for p in ${TARGET}; do 
    echo -e "\nMaking target $p" 
    make -j40 $p
done 
