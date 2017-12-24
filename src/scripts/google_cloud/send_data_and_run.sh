#!/bin/bash 

declare -a essential_files=(
    "config.xml"
    "kin"
    "models"
    "speaker_shader_assets"
    "launch_parallel_chunk.sh"
)

dest_dir="/home/jui-hsien/code/acoustics/work/demo/talk_fan_2"
create_dir="data"

echo -e "\n=== Sending files/dirs START ==="
for i in "${essential_files[@]}"; do
    echo $i;
done
echo -e "=== Sending files/dirs END ==="

for ii in `seq 7 26`; do 
    p="instance-${ii}"
    chunk_index=`expr ${ii} - 7`
    echo -e "\n\nSending data to ${p} and run wavesolver\n"
    gcloud compute ssh $p -- -t "mkdir -p ${dest_dir}; rm -r ${dest_dir}/${create_dir}; rm -r ${dest_dir}/models/*;" 
    for file in "${essential_files[@]}"; do
        gcloud compute scp --recurse ${file} $p:${dest_dir}/
    done
    gcloud compute ssh $p -- -t "
        cd ${dest_dir}; 
        mv ${create_dir} ${create_dir}_old;
        mkdir -p ${create_dir}; 
        tmux new -d \". ~/.zshrc;. ./launch_parallel_chunk.sh 32 config.xml ${chunk_index} run.log;\" "
done
         
