#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=run_all_pt2
#SBATCH --mem=100000
#SBATCH --time=1:20:00
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

# module load anaconda3/2023.09 # load modules... not needed since it's loaded by default
python scripts/DEBIAS-M_semi_leaky.py $1 $2 $3 &> ./output/$3/logs/run_DEBIAS-M_out_${1}.log &

wait

echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
