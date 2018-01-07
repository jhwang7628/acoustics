#!/bin/bash

declare -a essential_files=(
    "assets"
    "chunks.txt"
    "config.xml"
    "launch_parallel_chunk.sh"
)

dest_dir="/home/jui-hsien/code/acoustics/work/demo/character_poser_2"
create_dir="data"

echo -e "\n=== Sending files/dirs START ==="
for i in "${essential_files[@]}"; do
    echo $i;
done
echo -e "=== Sending files/dirs END ==="

for ii in `seq 8 24`; do
    p="instance-${ii}"
    chunk_index=`expr ${ii} - 7`
    echo -e "\n\nSending data to ${p} and run wavesolver\n"
    gcloud compute ssh $p -- -t "mkdir -p ${dest_dir}; rm -r ${dest_dir}/${create_dir}; rm -r ${dest_dir}/models/*;"
    for file in "${essential_files[@]}"; do
        gcloud compute scp --recurse ${file} $p:${dest_dir}/
    done
    gcloud compute ssh $p -- -t "
        cd ${dest_dir};
        pushd assets/Dilo01;
        if [ ! -d processed ]; then unzip processed.zip; fi;
        if [ ! -d sphere ]; then unzip sphere.zip; fi;
        popd;
        if [ -d ${create_dir} ]; then mv ${create_dir} ${create_dir}_old; fi;
        mkdir -p ${create_dir};
        tmux new -d \". ~/.zshrc;. ./launch_parallel_chunk.sh 32 config.xml ${chunk_index} run.log;\" "
done

