#!/bin/bash 
if [ $# -ne 3 ]; then
    echo "**Usage: ${0} <model_name> <tets_folder> <cp_destination>"
    exit
fi
MODEL_NAME=${1}
TETS_FOLDER=${2}
DEST=${3}
TARGET=${TETS_FOLDER}/${MODEL_NAME}.tet
if [ -f ${TARGET} ]; then
    echo "cp ${TARGET} ${DEST}" 
    cp ${TARGET} ${DEST}
else 
    echo " Cannot locate ${TARGET}" 
fi
