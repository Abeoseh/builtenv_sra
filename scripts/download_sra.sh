#!/bin/bash

#export PATH=$PATH:"/c/Program Files/sratoolkit.3.1.1-win64/bin"  # I added it as an environmental variable to my path 
# projects=("PRJEB14474")
projects=("PRJEB6292" "PRJEB33050" "PRJEB26708")

for project in ${projects[@]}; do
	echo "start" ${project}
	while read SRR  ; do
		## check if the fastq.gz file for that sample has been fetched
		if [ ! -f ../${project}/fastq/$SRR.fastq.gz ]; then 
			## if not then create the fastq.gz file
			echo $SRR
			prefetch $SRR
			fasterq-dump.exe $SRR -v --split-3 --outdir ../${project}/fastq
			gzip ../${project}/fastq/${SRR}*
			l="./$SRR"
			echo $l
			rm -rf $l
		fi
	done < ../${project}/SRR_Acc_List.txt
	echo "done" ${project}
	echo "_______________________________________"
done