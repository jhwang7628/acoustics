#!/bin/bash 

declare -a essential_files=(
    # "models"
    # "rush_working_man_mono.wav" 
    # "config_parallel.xml"
    "launch_parallel_chunk.sh"
)

dest_dir="/home/jui-hsien/code/acoustics/work/demo/musical_box_1"
create_dir="data_test"

echo -e "\n=== Sending files/dirs START ==="
for i in "${essential_files[@]}"; do
    echo $i;
done
echo -e "=== Sending files/dirs END ==="

for ii in `seq 8 26`; do 
    p="instance-${ii}"
    chunk_index=`expr ${ii} - 7`
    echo -e "\n\nSending data to ${p} and run wavesolver\n"
    gcloud compute ssh $p -- -t "mkdir -p ${dest_dir};" 
    for file in "${essential_files[@]}"; do
        gcloud compute scp --recurse ${file} $p:${dest_dir}/
    done
    gcloud compute ssh $p -- -t "
        cd ${dest_dir}; 
        mkdir -p ${create_dir}; 
        tmux new -s run -d \". ~/.zshrc;. ./launch_parallel_chunk.sh 32 config_parallel.xml ${chunk_index} run.log;\" "
done
         
