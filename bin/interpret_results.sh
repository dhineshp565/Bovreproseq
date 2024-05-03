#!/usr/bin/env bash

# This script checks for the presence of target names in the 15th column of the abricate output
# generates a csv file with four columns the results (positive or negative)
# Also provides warning if the internal amplification control is negative
# $1 = SampleName, $2 = List of target names

# start by creating the header of the csv file
echo -e "SAMPLEID,PATHOGEN,RESULT,STATUS" > "$1_results.csv"

# loop through the list of target names
while read lines; do
	target=$(echo "${lines}" | cut -f 1) # extract the target name
    amplicon=$(cut -f 15 "$1_abricate.csv"| grep -F "${target}") # check for the presence of target name in the 15th column of the abricate output
    if [ "${amplicon}" ]; then
        echo -e "$1,${target},Positive,OK" >> "$1_results.csv" # if the target name is present
    else
        echo -e "$1,${target},Negative,OK" >> "$1_results.csv"  # if the target name is not present
    fi
done < $2

# check for the presence of Internal_amplification_control,Negative in the csv file
if grep -Fq "Internal_amplification_control,Negative" "$1_results.csv"; then
    sed -i 's/OK/Check for PCR inhibition/g' "$1_results.csv" # if the Internal_amplification_control is negative then change status
fi
