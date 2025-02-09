#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --time=200:00:00  
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=20GB
#SBATCH -o ../logs/qiime2/%a_paired.out # STDOUT

## this file takes one positional argument: $1 which is the folder name containing the fastq folder and manifest file, it is assumed that the folder is the SRA project name

module load qiime2/2024.2 

## remove dada2 folder if it exists
rm -rf /users/aflemis1/sink/$1/qiime2_output/dada2


qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /users/aflemis1/sink/$1/manifest.tsv \
  --output-path /users/aflemis1/sink/$1/qiime2_output/paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 

qiime dada2 denoise-paired \
	--i-demultiplexed-seqs /users/aflemis1/sink/$1/qiime2_output/paired-end-demux.qza \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-trunc-len-f 0 \
	--p-trunc-len-r 0 \
	--p-chimera-method pooled \
	--p-min-fold-parent-over-abundance 100 \
	--p-n-threads 30 \
	--output-dir /users/aflemis1/sink/$1/qiime2_output/dada2


qiime feature-classifier classify-sklearn \
      --i-classifier /users/aflemis1/Taxonomy_qiime2/silva-138-99-nb-classifier.qza \
      --p-n-jobs 1 \
      --i-reads /users/aflemis1/sink/$1/qiime2_output/dada2/representative_sequences.qza \
      --o-classification /users/aflemis1/sink/$1/qiime2_output/dada2/taxonomy.qza

qiime tools export \
      --input-path  /users/aflemis1/sink/$1/qiime2_output/dada2/representative_sequences.qza \
      --output-path  /users/aflemis1/sink/$1/qiime2_output/dada2/representative_sequences

qiime tools export \
      --input-path  /users/aflemis1/sink/$1/qiime2_output/dada2/taxonomy.qza \
      --output-path  /users/aflemis1/sink/$1/qiime2_output/dada2/taxonomy

qiime tools export \
      --input-path /users/aflemis1/sink/$1/qiime2_output/dada2/table.qza \
      --output-path  /users/aflemis1/sink/$1/qiime2_output/dada2/table
      
qiime tools export \
      --input-path /users/aflemis1/sink/$1/qiime2_output/dada2/denoising_stats.qza  \
      --output-path  /users/aflemis1/sink/$1/qiime2_output/dada2/denoising_stats

## merge taxonomy with the feature table
biom add-metadata \
 --input-fp /users/aflemis1/sink/$1/qiime2_output/dada2/table/feature-table.biom \
 --observation-metadata-fp /users/aflemis1/sink/$1/qiime2_output/dada2/taxonomy/taxonomy.tsv \
 --output-fp /users/aflemis1/sink/$1/qiime2_output/dada2/taxonomy/biom-with-taxonomy.biom  \
 --sc-separated taxonomy --observation-header OTUID,taxonomy 




biom convert -i /users/aflemis1/sink/$1/qiime2_output/dada2/taxonomy/biom-with-taxonomy.biom -o ../csv_files/combine/$1.txt --to-tsv --header-key taxonomy

 
 
  




































