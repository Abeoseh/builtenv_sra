#!/bin/bash 
#SBATCH --partition=Orion
#SBATCH --job-name=combine_pt2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=100000
echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo "" cd $SLURM_SUBMIT_DIR


if [ ! -d "./logs/fix_taxonomy" ]
then
	mkdir ./logs/fix_taxonomy
fi

module load anaconda3/2023.09 # load modules


python ./scripts/fix_taxonomy.py &> logs/fix_taxonomy/fix_taxonomy.log


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
