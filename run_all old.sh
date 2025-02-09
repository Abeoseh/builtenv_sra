#!/bin/bash

## to run this file: 
# run_all.sh input_folder_name output_folder_name amount_of_studies
# run_all.sh associated associated 4

# If the folder already exists allow user to decide to overwrite
while [ -d ./output/$2 ]; do
	read -p "This folder exists. Overwrite? Y (yes) or N (no) " yn
		case $yn in
			[Yy]* ) rm -r ./output/$1; break;;
			[Nn]* ) exit;;
			* ) echo "Please answer Y (yes) or N (no)";;
		esac
done


mkdir ./output/$2
mkdir ./output/$2/logs
mkdir ./output/$2/AUCs
mkdir ./output/$2/DEBIAS-M_runs
mkdir ./output/$2/ROC_histograms


search_dir=./output/$2/DEBIAS-M_runs
amount=$3


for i in $( eval echo {1..$amount} )
do
	sbatch run_all_pt1.sbatch $i 100 $1 $2
	sbatch run_all_pt2.sbatch $i $1 $2
	echo $i "of $amount done"
done



count_files() {
	file_count=$(find "$search_dir/"*debiased* | wc -l)
	echo $file_count
}

while [ $(count_files) -lt $amount ]; do
	echo "Waiting for $amount files"
	echo $file_count
	sleep 150 # waiting 2 minute 30s before checking again
done

echo "$amount files have been detected"


# entry is each DEBIASed file
for entry in "$search_dir/"*debiased*
do
	echo "$entry"
	# extract study IDs
	Study_ID="${entry##*/}"  # remove everything before the last /
	Study_ID="${Study_ID##*_}"  # Remove everything before the last _
	Study_ID="${Study_ID%.*}"    # Remove everything after the last .

	# Untested ALT
	# echo "$file" | grep -oP '([^_]+)(?=\.\w+$)' 

	sbatch run_all_pt3.sbatch $entry 100 $Study_ID $1

	echo "$Study_ID"
done

echo "done with all"



