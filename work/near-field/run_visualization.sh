#!/bin/bash

if [ $# -ne 2 ]; then 
    echo "**Usage: ${0} <work_dir> <max_number_files>"
    exit
fi

WORKDIR=${1}
XML_FILE=$(ls ${1}/*.xml)
/home/jui-hsien/code/acoustics/build/bin/convert-dat-vtk ${XML_FILE} ${1} ${2}
