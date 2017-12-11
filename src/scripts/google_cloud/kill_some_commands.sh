#!/bin/bash 

if [ $# -ne 1 ]; then 
    echo "**Usage: ${0} <kill_command_name_prefix>" 
else 
    kill_command_name=${1}
    for ii in `seq 7 26`; do
        p="instance-${ii}" 
        gcloud compute scp kill_process.sh ${p}:~
        gcloud compute ssh ${p} -- -t ". ./kill_process.sh ${kill_command_name};"
    done
fi


