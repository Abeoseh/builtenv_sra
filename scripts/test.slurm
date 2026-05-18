#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --time=200:00:00  
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=20GB
#SBATCH -o ../logs/qiime2/%j_paired.out # STDOUT

## this file takes one positional argument: $1 which is the folder name containing the fastq folder and manifest file, it is assumed that the folder is the SRA project name

module load qiime2/2024.2 
cd ../PRJEB11111/qiime2_output
pwd 

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv

echo "done"

