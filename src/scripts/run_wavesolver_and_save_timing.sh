#!/bin/zsh

if [ $# -ne 3 ]; then 
    echo "**Usage: ${0} <num_threads> <config> <log_file>"
    exit 1
fi
echo -e "\nRunning wavesolver: `pwd`"
echo -e " OMP_NUM_THREADS: ${1}" 
echo -e " config         : ${2}" 
echo -e " log file       : ${3}"
{time OMP_NUM_THREADS=${1} fdtd-acoustic-simulator-viewer --config ${2} --nogui } >${3} 2>&1
