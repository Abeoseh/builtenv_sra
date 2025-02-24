**mk_manifest.sh $1 $2**
```
How to run:
./scripts/mk_manifest.sh fastq_folder output_folder

ex:
./scripts/mk_manifest.sh ../PRJEB3232/fastq ../PRJEB3232
assumes the directory structure ./SRA_project/fastq/*.fastq
can be used for single or paired runs 
```

**qiime2_single.slurm**
```
single means only forward or reverse
ensure folders: qiime2_output and ./

How to run:
sbatch scripts/qiime2_single.slurm folder_containing_fastq_folder

ex:
sbatch scripts/qiime2_single.slurm PRJEB3232
```

**qiime2_paired.slurm**
```
paired means forward and reverse

How to run:
sbatch scripts/qiime2_paired.slurm folder_containing fastq folder

ex:
sbatch scripts/qiime2_paired.slurm PRJNA834026
```


**./scripts/combine_pt3.R**
```
$1 output folder name
$2 first phenotype
$3 second phenotype
$4 (""/T) keep NA as NA
./scripts/combine_pt3.R $1 $2 $3 $4

ex:
sbatch combine_pt3.sbatch "associated" "skin associated" "floor associated" # results in NA values being coded as 0
sbatch combine_pt3.sbatch "associated_na" "skin associated" "floor associated" "T" # results in NA staying as NA
```
**./scripts/PERMANOVA_PCOA.R**
```
Makes PERMANOVA and POCA plots
$1 full input filepath including file name
$2 full output filepath including filename without .png, but no error will occur if .png is given.
$3 is pheno1
$4 is pheno2
./scripts/PERMANOVA_PCOA.R $1 $2 $3 $4

ex: 
sbatch PERMANOVA_PCOA.sbatch ./csv_files/vetted_ontology/lognorm_data.csv ./output/PCOA/vetted
sbatch PERMANOVA_PCOA.sbatch ./csv_files/vetted_ontology/lognorm_data.csv ./output/PCOA/vetted.png # the .png will be stripped of
```

**./scripts/pval-pval_plot.R**
```
Makes pval-pval plots
sbatch pval-pval_plot.sbatch $1 $2
$1 input folder with lognorm_data.csv inside
$2 output folder located within ./output

ex:
sbatch pval-pval_plot.sbatch sink_nowater sink_nowater
```
