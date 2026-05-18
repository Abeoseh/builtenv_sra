All initial fastq downloading was done locally. Afterwards, files were processed on the cluster.

skin vs floor skin is 1 and floor is 0<br>
skin vs skin associated: skin is 1 and skin associated is 0<br>
skin associated vs floor: skin associated is 1 and floor is 0<br>

*_longitudinal includes all timepoints, including ones where humans didn't inhabit that building (ex pre-opening for the hospital)<br>
while the *_timepoints don't include pre-inhabiting timepoints


**download_sra.sh**
```
`download_sra.sh` takes a SRR_Acc_List.txt and downloads all the fastq files for the given accessions. It automatically
loop over the project IDs provided in as an array called "projects" within `download_sra.sh`.

How to run:
- Change the projects in the array named "projects" to your desired project(s). 
        - The script assumes you named your folders with the project name. All of this can be changed.
- Open command prompt and navigate to the scripts folder.
- Run the following command: `download_sra.sh > sra.log 2>&1` 
- Output: a file in the scripts folder called sra.log which contains std error and std Output.
```

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

**fix_taxonomy.py**
```
assumes the file name (ex PRJEB11111.txt) is the same as the folder name with the qiime2 manifest.tsv

How to run:
scripts/fix_taxonomy.py
output: single_*study name*.csv for all the studies

```

**./scripts/combine_pt2.R**
```
combines all the files from fix_taxonomy.py
output combine_otus.csv

ex:
sbatch combine_pt2.sbatch
```

**./scripts/combine_pt3.R**
```
$1 output folder name
$2 first phenotype
$3 second phenotype
$4 ""/T keep NA as NA
$5 ""/T add a timepoint column
./scripts/combine_pt3.R $1 $2 $3 $4 $5

ex:
sbatch combine_pt3.sbatch "associated" "skin associated" "floor associated" # results in NA values being coded as 0
sbatch combine_pt3.sbatch "associated_na" "skin associated" "floor associated" "T" # results in NA staying as NA
sbatch combine_pt3.sbatch skin_floor_timepoint skin floor "" T # results as NA values being coded as 0 ($4="") and a timepoint column added ($5=T)
```

**./scripts/PERMANOVA.R**
```
Makes PERMANOVA plots for all samples, hand, and environmental samples. 
$1 full input filepath including file name
$2 full output filepath including filename without .png, but no error will occur if .png is given.
$3 is pheno1 (hand)
$4 is pheno2 (environmental)
./scripts/PERMANOVA_PCOA.R $1 $2 $3 $4

ex: 
sbatch PERMANOVA_PCOA.sbatch ./csv_files/vetted_ontology/lognorm_data.csv ./output/PCOA/vetted
sbatch PERMANOVA_PCOA.sbatch ./csv_files/vetted_ontology/lognorm_data.csv ./output/PCOA/vetted.png # the .png will be stripped of
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

**./run_all.sh**
```
$1 input folder (ensure it's located within ./csv_files)
$2 output folder (ensure folder ./output is created)
$3 amount of studies
./run_all.sh $1 $2 $3

ex:
./run_all.sh sink_nonsink sink_nonsink 4
```

**timepoint PCoAs**
makes:
- bar chart of R^2 for the PERMANOVA for sample type, time points, and ontology (original sample names before we merged them into hand, hand associated, and floor)
  - NOTE: the ontology column is automatically created.
- PCoAs for each time point colored by sample type
- PCOA for each study colored by all time points 
- PCOA for all studies together colored by sample name
- PCOA between 0013/07 and other timepoints colored by phenotype (sample name) and timepoint

``` 
$1 input lognorm file
$2 output path with file name (no extension)
$3 pheno1 
$4 pheno2 
$5 Whether to do the PERMANOVAs again. There isn't a need to recalcualte them if they were already done. (/T)

ex:
sbatch timepoint_pcoa.sbatch ./csv_files/skinassoc_timepoint/lognorm_data_all.csv \
./output/skinassoc_timepoint/PCOA/skinassoc hand "hand associated" # Doesn't recalculate the PERMANOVAs

```


Processing Notes:
- PRJEB3232 and PRJEB3250: only have one read per spot.
- PRJNA834026 had no taxonomy assigned and was dropped.


 ### provide accessions:
    ### automatically makes folders then downloads data (run table... then renames to project_run_table and accession list)
