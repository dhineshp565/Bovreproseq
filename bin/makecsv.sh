#!/bin/env bash


# get list of directories and replaced slash
ls -1 $1 > samples.csv

# get path of each folder
realpath $1/* > paths.csv

# concatenate samplenames and path with comma as delimiter
paste samples.csv paths.csv > samplelist.csv
sed -i 's/	/,/g' samplelist.csv
	
# add headers to the csv file
sed -i '1i SampleName,SamplePath' samplelist.csv