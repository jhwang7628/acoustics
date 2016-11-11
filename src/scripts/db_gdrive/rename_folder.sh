#!/bin/bash 
FOLDER=downloaded
REDUNDENT_PREFIX=p_
for p in ${FOLDER}/*; do 
    FILENAME=$(basename "$p")
    MODEL_NAME="${FILENAME##*${REDUNDENT_PREFIX}}"
    if [ ! -d ${MODEL_NAME} ]; then
        echo "cp -r ${p} ${MODEL_NAME}"
        cp -r ${p} ${MODEL_NAME}
    fi
done
