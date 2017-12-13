#!/bin/bash 

BUILD_DIR="/home/jui-hsien/code/acoustics/build_release" 
N_CORES=32

for ii in `seq 7 26`; do 
    p="instance-$ii"
    echo -e "\n ============== $p ===============" 
    echo "Launching make on remote node and detach"
    gcloud compute ssh $p -- -t "
        cd ${BUILD_DIR}; 
        git log | head -6
    "
done
