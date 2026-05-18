#!/bin/bash

## to run this file: 
# run_all.sh input_folder_name output_folder_name amount_of_studies
# run_all.sh associated associated 4

# If the folder already exists allow user to decide to overwrite

remove_dirs=("./output/$2/CSVs" "./output/$2/plots")
while [ -d ./output/$2/CSVs ]; do
	read -p "This script has been run before. Remove CSVs and plots? Y (yes) or N (no) " yn
		case $yn in
			[Yy]* ) 
				for dir in "${remove_dirs[@]}"; do
					rm -r $dir;
				done
			break;;

			[Nn]* ) exit;;

			* ) echo "Please answer Y (yes) or N (no)";;
		esac
done


# remove_dirs=("./output/$2/DEBIAS-M_runs")
if [ -d "./output/$2/DEBIAS-M_runs" ]; then
	read -p "Remove DEBIAS-M files? Y (yes) or N (no) " yn
		case $yn in
			[Yy]* ) 
				rm -r "./output/$2/DEBIAS-M_runs"
				;;
			[Nn]* ) echo "no DEBIAS-M files removed"
				;;
			* ) 
				echo "Please answer Y (yes) or N (no)"
				;;
		esac
fi



mkdir ./output/$2
mkdir ./output/$2/logs
mkdir ./output/$2/CSVs
mkdir ./output/$2/DEBIAS-M_runs
mkdir ./output/$2/plots


search_dir=./output/$2/DEBIAS-M_runs
amount=$3

batches=100
## submit batches of 100, 10 times for a total of 1000 permutations

## run pre-Random forest

for i in $( eval echo {1..$amount} )
do
	for batch in {1..10}
	do
		sbatch run_all_pt1.sbatch $i 100 $1 $2 $batch
		echo "batch: $batch for study $i done"
	done

	echo $i "of $amount done pre-random forest"
	echo "_______________________________________________"
done



## run DEBIAS-M 
for i in $( eval echo {1..$amount} )
do
	sbatch run_all_pt2.sbatch $i $1 $2
	echo $i "of $amount submitted to DEBIAS-M"
done

echo "_______________________________________________"

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

#	# Untested ALT
#	# echo "$file" | grep -oP '([^_]+)(?=\.\w+$)' 
	for batch in {1..10}
	do
		sbatch run_all_pt3.sbatch $entry 100 $Study_ID $2 $batch

	done

	echo "$Study_ID done post random forest"
	echo "_______________________________________________"
done

echo "done with all"



