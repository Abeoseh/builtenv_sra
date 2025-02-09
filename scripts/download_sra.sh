#!/bin/bash

export PATH=$PATH:"/c/Program Files/sratoolkit.3.1.1-win64/bin"
projects=("PRJNA878661" "PRJNA834026" "PRJEB3250" "PRJEB3232")
for project in ${projects[@]}; do
	echo "start" ${project}
	while read SRR  ; do
        echo $SRR
        fasterq-dump.exe $SRR -v --split-3 --outdir ../${project}/fastq
        gzip ../${project}/fastq/${SRR}*
	done < ../${project}/SRR_Acc_List.txt
	echo "done" ${project}
	echo "_______________________________________"
done