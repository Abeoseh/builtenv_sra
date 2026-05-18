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
					printf "not removing plot and CSV directories"
				done
			break;;

			[Nn]* ) exit;;

			* ) echo "Please answer Y (yes) or N (no)";;
		esac
done




mkdir ./output/$2
mkdir ./output/$2/logs
mkdir ./output/$2/CSVs
mkdir ./output/$2/plots


search_dir=./output/$2/DEBIAS-M_runs
amount=$3

batches=100
## submit batches of 100, 10 times for a total of 1000 permutations

## run pre-Random forest

for i in $( eval echo {1..$amount} )
do
	for batch in {0..3}
	do
		sbatch run_all_timepoint_pt1.sbatch $i 3 $1 $2 $batch
		echo "batch: $batch for study $i done"
	done

	echo $i "of $amount done pre-random forest"
	echo "_______________________________________________"
done





