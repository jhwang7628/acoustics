#!/bin/bash 

ANALYTICAL_FOLDER=analytical
SOUND_FOLDER=sound

mkdir ${ANALYTICAL_FOLDER}
for p in test_*analytical*; do 
    echo $p; 
    mv $p ${ANALYTICAL_FOLDER} 
done

mkdir ${SOUND_FOLDER}
mv *.wav ${SOUND_FOLDER}
