#!/bin/bash

BUILD_DIR="/home/jui-hsien/code/acoustics/build_release"
N_CORES=16

for ii in `seq 7 24`; do
    p="instance-$ii"
    echo -e "\n ============== $p ==============="
    echo "Launching make on remote node and detach"
    gcloud compute ssh $p -- -t "
        cd ${BUILD_DIR};
        git pull
        tmux new -s make -d \"
            . ~/.zshrc;
            cmake ..;
            make fdtd-acoustic-simulator-viewer -j ${N_CORES} > make.log 2>&1;\"
    "
done
