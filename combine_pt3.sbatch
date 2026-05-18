#!/bin/bash 
#SBATCH --partition=Orion
#SBATCH --job-name=combine_pt3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100000
echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo "" cd $SLURM_SUBMIT_DIR

if [ ! -d ./csv_files/$1 ]; then
	mkdir ./csv_files/$1
fi

if [ "$4" == "" ]; then logs=logs/combine/combine_pt3_$1.log
else logs=logs/combine/combine_pt3_$1_na.log; fi

module load R

echo "$ 1: $1"
echo "$ 2: $2"
echo "$ 3: $3"
echo "$ 4: $4"
echo "$ 5: $5"

Rscript ./scripts/combine_pt3.R "$1" "$2" "$3" "$4" "$5" &> $logs


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
