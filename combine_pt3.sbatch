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

module load R

Rscript ./scripts/combine_pt3.R "$1" "$2" "$3" $4 &> logs/combine/combine_pt3_$1.log


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
