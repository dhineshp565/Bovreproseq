#!/usr/bin/env bash
# This script differentiates Campylobacter fetus subspecies into .fetus or .venerealis based on ST typing (field #3,row #2 in MLST ouput)
# appends the subspecies resultto the the last column of the results file
#if no ST type is available NA is added to the last column

# Campylbacter fetus venerealis typing
ST="$(tail -n +2 $1 | cut -f 3)"
if [[  $ST =~ ^(4|7|12) ]];then
	subspecies="Campylobacter fetus. venerealis" 
# Campylbacter fetus fetus typing	
elif [[ $ST =~ ^(1|2|3|5|6|8|9|10|11|13|14)$ ]];then
	subspecies="Campylobacter fetus. fetus"
# handle No ST typing		
else 
	subspecies="NA"		
fi
# add results
paste "$1" <(printf "ORGANISM\n%s" "$subspecies") > $(basename $1 .csv)_results.csv