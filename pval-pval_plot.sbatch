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


if [ ! -d "./output/$2/pval_pval" ]; then
	mkdir ./output/$2/pval_v_pval

fi

if [ ! -d "./output/$2/logs" ]; then
	mkdir ./output/$2/logs

fi


Rscript ./scripts/pval-pval_plot.R $1 $2 &> ./output/$2/logs/pval-pval_plot.log


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"


