#!/bin/bash

#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=run_all_pt3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100000
#SBATCH --time=48:20:00

echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
cd $SLURM_SUBMIT_DIR


module load R # load modules

Rscript scripts/random_forest_2.R $1 $2 $3 $4 &> ./output/$4/logs/run_rf_post_$3.log



echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
