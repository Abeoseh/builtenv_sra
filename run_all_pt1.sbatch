#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=run_all_pt1
#SBATCH --mem=100000
#SBATCH --time=28:20:00
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --nodes 2
echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
cd $SLURM_SUBMIT_DIR

module load R # load modules

Rscript scripts/random_forest_1.R $1 $2 $3 $4 &> ./output/$4/logs/run_rf_pre_$1.log 

echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
