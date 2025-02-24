# Purpose:

This project is an integrated analysis between hand and floor samples to gleam insights about the microbiome which distinguishes them. All fastq files are publically available on the SRA.

All initial fastq downloading was done locally. Afterwards, files were processed on the cluster.


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

**./run_all.sh**
```
$1 input folder (ensure it's located within ./csv_files)
$2 output folder (ensure folder ./output is created)
$3 amount of studies
./run_all.sh $1 $2 $3
ex:
./run_all.sh sink_nonsink sink_nonsink 4
```

Processing Notes:
- PRJEB3232 and PRJEB3250: only have one read per spot.
- PRJNA834026 had no taxonomy assigned and was dropped.


 ### provide accessions:
    ### automatically makes folders then downloads data (run table... then renames to project_run_table and accession list)
