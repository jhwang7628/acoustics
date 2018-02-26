#!/bin/bash

declare -a essential_files=(
    "assets"
    "imp_fullspace_sim-014.06.51.tar"
    "fdtd-transfer-cloud"
)

dest_dir="/home/jui-hsien/code/shells/work/cymbal_drum_stick_2"
create_dir="data_14"

echo -e "\n=== Sending files/dirs START ==="
for i in "${essential_files[@]}"; do
    echo $i;
done
echo -e "=== Sending files/dirs END ==="

for ii in `seq 11 19`; do
    p="instance-${ii}"
    chunk_index=`expr ${ii} - 10`
    echo -e "\n\nSending data to ${p} and run wavesolver\n"
    gcloud compute ssh $p -- -t "mkdir -p ${dest_dir};"
    for file in "${essential_files[@]}"; do
        gcloud compute scp --recurse ${file} $p:${dest_dir}/
    done
    gcloud compute ssh $p -- -t "
        cd ${dest_dir};
        tar -xvf imp_fullspace_sim-014.06.51.tar;
        cd fdtd-transfer-cloud;
        if [ -d ${create_dir} ]; then rm -r ${create_dir}_old; mv ${create_dir} ${create_dir}_old; fi;
        mkdir -p ${create_dir};
        tmux new -d \". ~/.zshrc;. ./launch_parallel_chunk.sh 32 config.xml ${chunk_index} run_${chunk_index}.log;\" "
done

