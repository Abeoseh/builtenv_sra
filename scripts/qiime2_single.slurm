#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --time=200:00:00  
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=20GB
#SBATCH -o ../logs/qiime2/%j_single.out # STDOUT

echo "Start Time : $(date)"
echo "folder: ${1}"

## this file takes one positional argument: $1 which is the folder name containing the fastq folder and manifest file, it is assumed that the folder is the SRA project name
## folder qiime2_output needs to be made within $1

module load qiime2/2024.2 

## remove dada2 folder if it exists
rm -rf /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2


qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /users/aflemis1/builtenv_sra/$1/manifest.tsv \
  --output-path /users/aflemis1/builtenv_sra/$1/qiime2_output/single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 



qiime dada2 denoise-single \
	--p-trim-left 0 \
	--p-trunc-len 0 \
	--i-demultiplexed-seqs /users/aflemis1/builtenv_sra/$1/qiime2_output/single-end-demux.qza \
	--p-chimera-method pooled \
	--p-min-fold-parent-over-abundance 100 \
	--p-n-threads 30 \
	--output-dir /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2


qiime feature-classifier classify-sklearn \
      --i-classifier /users/aflemis1/Taxonomy_qiime2/silva-138-99-nb-classifier.qza \
      --p-n-jobs 1 \
      --i-reads /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/representative_sequences.qza \
      --o-classification /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/taxonomy.qza

qiime tools export \
      --input-path  /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/representative_sequences.qza \
      --output-path  /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/representative_sequences

qiime tools export \
      --input-path  /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/taxonomy.qza \
      --output-path  /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/taxonomy

qiime tools export \
      --input-path /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/table.qza \
      --output-path  /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/table
      
qiime tools export \
      --input-path /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/denoising_stats.qza  \
      --output-path  /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/denoising_stats

## merge taxonomy with the feature table
biom add-metadata \
 --input-fp /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/table/feature-table.biom \
 --observation-metadata-fp /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/taxonomy/taxonomy.tsv \
 --output-fp /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/taxonomy/biom-with-taxonomy.biom  \
 --sc-separated taxonomy --observation-header OTUID,taxonomy 




biom convert -i /users/aflemis1/builtenv_sra/$1/qiime2_output/dada2/taxonomy/biom-with-taxonomy.biom -o ../csv_files/combine/$1.txt --to-tsv --header-key taxonomy

 
 
  
echo "script complete, for assigned taxnomy file check ../csv_files/combine/ for ${1}.txt"
echo "End Time : $(date)"