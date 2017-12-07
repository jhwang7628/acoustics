#!/bin/bash

cd ${HOME}/code/acoustics/build_release
rm -r *
cmake ..
make fdtd-acoustic-simulator-viewer -j16

cd ${HOME}/code/shells
. ./.setup.sh
make clean
make 
