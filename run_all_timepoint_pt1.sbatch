#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=run_all_pt1
#SBATCH --mem=100000
#SBATCH --time=96:20:00
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

Rscript scripts/timepoint_random_forest.R $1 $2 $3 $4 $5 &> ./output/$4/logs/run_rf_pre_$1_batch_$5.log 

echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
