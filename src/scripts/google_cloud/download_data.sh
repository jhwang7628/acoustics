#!/bin/bash

declare -a essential_files=(
    "/home/jui-hsien/code/acoustics/work/demo/character_poser_2/data/*_start_time*"
    "/home/jui-hsien/code/acoustics/work/demo/character_poser_2/data/*_metadata*"
    "/home/jui-hsien/code/acoustics/work/demo/character_poser_2/data/*_listening_position.dat"
    "/home/jui-hsien/code/acoustics/work/demo/character_poser_2/data/*_all_audio.dat"
)

dest_dir="data"
mkdir -p ${dest_dir}

for ii in `seq 7 24`; do
    p="instance-${ii}"
    for file in "${essential_files[@]}"; do
        echo "gcloud compute scp --recurse $p:${file} ${dest_dir}/"
        gcloud compute scp --recurse $p:${file} ${dest_dir}/
    done
done

