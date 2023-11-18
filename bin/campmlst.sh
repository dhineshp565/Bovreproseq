#!/bin/env bash
# This script differentiates Campylobacter fetus subspecies into .fetus or .venerealis based on ST typing (field #3,row #2 in MLST ouput)
# appends the subspecies resultto the the last column of the results file
#if no ST type is available NA is added to the last column

# Campylbacter fetus venerealis typing
ST="$(tail -n +2 $1 | cut -f 3)"
if [  $ST -eq "4" -o $ST -eq "7" -o $ST -eq "12" ];then
	echo "ORGANISM" > temp.csv
	echo "Campylobacter fetus. venerealis" >> temp.csv
	paste $1 temp.csv > $(basename $1 .csv)_results.csv
# Campylbacter fetus fetus typing	
elif [ $ST -eq "1" -o $ST -eq "2" -o $ST -eq "3" -o $ST -eq "5" -o $ST -eq "6" -o $ST -eq "8" -o $ST -eq "9" -o $ST -eq "10" -o $ST -eq "11" -o $ST -eq "13" -o $ST -eq "14" ];then
	echo "ORGANISM" > temp.csv
	echo "Campylobacter fetus. fetus" >> temp.csv
	paste $1 temp.csv > $(basename $1 .csv)_results.csv
# handle No ST typing		
else 
	echo "ORGANISM" > temp.csv
	echo "NA" >> temp.csv
	paste $1 temp.csv > $(basename $1 .csv)_results.csv
		
fi