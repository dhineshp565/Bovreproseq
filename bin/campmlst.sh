#!/bin/env bash
shopt -s extglob
	ST="$(tail -n +2 $1 | cut -f 3)"
	if [  $ST -eq "4" -o $ST -eq "7" -o $ST -eq "12" ];then
		echo "ORGANISM" > temp.csv
		echo "Campylobacter fetus. venerealis" >> temp.csv
		paste $1 temp.csv > $(basename $1 .csv)_results.csv
	
	elif [ $ST -eq "1" -o $ST -eq "2" -o $ST -eq "3" -o $ST -eq "5" -o $ST -eq "6" -o $ST -eq "8" -o $ST -eq "9" -o $ST -eq "10" -o $ST -eq "11" -o $ST -eq "13" -o $ST -eq "14" ];then
		echo "ORGANISM" > temp.csv
		echo "Campylobacter fetus. fetus" >> temp.csv
		paste $1 temp.csv > $(basename $1 .csv)_results.csv
		
	else 
		
		echo "ORGANISM" > temp.csv
		echo "NA" >> temp.csv
		paste $1 temp.csv > $(basename $1 .csv)_results.csv
		
	fi