#!/bin/bash

projects=("PRJEB14474")

# makes SRR_Acc_List from existing SRR_Acc_List without entries that have already been read


while read line; do
	if [ ! -f ../PRJEB14474/fastq/$line.fastq.gz ]; then
		echo $line >> ../PRJEB14474/SRR_Acc_List.txt #this works
	fi
done < ../PRJEB14474/SRR_Acc_List_og.txt

