#!/bin/bash 
#SBATCH --partition=Orion
#SBATCH --job-name=pval-pval
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

# $1 is full input path including file name
# $2 is full output path including file name

Rscript ./scripts/PCOA.R $1 $2 &> ./logs/PCOA/PCOA.log


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"


