#!/bin/bash 

if [ $# -ne 2 ]; then 
    echo "**Usage: ${0} <source_dir> <target_dir>"
    exit
fi

SRC_DIR=${1}
TAR_DIR=${2}
FILE_LIST="*.wav *_listening_position.dat *_input_control_file.xml *_solver_settings.txt"
for f in ${FILE_LIST}; do
    echo "Moving file: ${f}"
    cp ${SRC_DIR}/${f} ${TAR_DIR}
done
