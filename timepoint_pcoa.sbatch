#!/bin/bash 
#SBATCH --partition=Orion
#SBATCH --job-name=Permanova_PCOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=100000
echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo "" cd $SLURM_SUBMIT_DIR

if [ ! -d "./logs/PERMANOVA_PCOA" ]
then
	mkdir ./logs/PERMANOVA_PCOA
fi



module load R

# $1 is full input path including file name
# $2 is full output path including file name
# $3 is pheno1
# $4 is pheno2
# $5 when supplied makes the permanovas

Rscript ./scripts/timepoint_pcoa.R $1 $2 "$3" "$4" $5 &> ./logs/PERMANOVA_PCOA/PERMANOVA_PCOA_timepoint_"$3"_"$4".log


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"


