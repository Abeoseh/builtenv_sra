#!/bin/bash

## to run this file: 
# run_all.sh input_folder_name output_folder_name amount_of_studies
# run_all.sh associated associated 4

echo "This script only runs random forest before DEBIAS-M"

# If the folder already exists allow user to decide to overwrite
remove_dirs=("./output/$2/AUCs" "./output/$2/DEBIAS-M_runs" "./output/$2/ROC_histograms")
while [ -d ./output/$2/AUCs ]; do
	read -p "This script has been run before. Remove output? Y (yes) or N (no) " yn
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



mkdir ./output/$2
mkdir ./output/$2/logs
mkdir ./output/$2/AUCs
mkdir ./output/$2/DEBIAS-M_runs
mkdir ./output/$2/ROC_histograms

amount=$3


for i in $( eval echo {1..$amount} )
do
	sbatch run_all_pt1.sbatch $i 100 $1 $2
	echo $i "of $amount done"
done


echo "done with all"



