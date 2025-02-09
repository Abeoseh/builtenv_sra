#!/bin/bash

# print the folder name
echo "positional argument 1 (input directory with fastq files): $1"
echo "positional argument 2 (output directory):" $2

# remove the manifest if it exists
if [ -f $2/manifest.tsv ] ; then
	rm $2/manifest.tsv
fi


i=1
for fastq in "$1"/*; do
	
	# check if forward and reverse reads are present
	if [[ "$fastq" == *"_"* ]]; then
		forward="True"
		# add file header
		if [[ $i == 1 ]]; then
			printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > $2/manifest.tsv

		fi
		# add sample ID and forward read
		if [[ "$fastq" == *"_1"* ]]; then
			sample=""
			sample="\nsample-$i\t$(realpath $fastq)\t" 
			i=$(($i+1))

		## add reverse read
		else
			
			sample="$sample$(realpath $fastq)"
			# echo $sample
			printf "$sample"  >> $2/manifest.tsv
			
			# i=$(($i+1)) this skips the first file for some reason

		fi
	fi
	# if only forward reads
	if [[ $forward == "" ]]; then
		
		# add file header
		if [[ $i == "1" ]] ; then
			printf "sample-id\tabsolute-filepath\t" > $2/manifest.tsv

		fi
		## add sample IDs and forward reads
		printf "\nsample-$i\t$(realpath $fastq)" >> $2/manifest.tsv
		# echo $(realpath $fastq)
		i=$(($i+1))


	fi
	
	# echo $i
	# echo "$fastq"
done


echo "finished"

