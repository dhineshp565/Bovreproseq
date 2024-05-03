#!/usr/bin/env bash

# This script makes a csv file with sample name and sample path with headers (SampleName,SamplePath).
# $1 = input path of fastq directories

# add headers to the csv file
echo "SampleName,SamplePath" >> samplelist.csv
input_dir="$1"
# iterate over a given input path
for dir in "$input_dir"/*; do
# check only for directories
    if [ -d "$dir" ];then
        samplename=$(basename "$dir")
        path=$(realpath "$dir")
        echo "${samplename},${path}" >> samplelist.csv
    fi
done