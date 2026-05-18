#!/bin/bash

#srun -n 1 --time=3:00:00 --partition=Nebula --nodes=1 --mem-per-cpu=8GB --pty /bin/bash

module load R

echo "beginning"

sbatch timepoint_pcoa.sbatch ./csv_files/associated_timepoint/lognorm_data_all.csv ./output/associated_timepoint/PCOA/associated \
"hand associated" floor


#Rscript ./scripts/timepoint_pcoa.R ./csv_files/associated_timepoint/lognorm_data_all.csv ./output/associated_timepoint/PCOA/associated \
#"hand associated" floor &> "./logs/PERMANOVA_PCOA/PERMANOVA_PCOA_timepoint_hand associated_floor associated.log"


## -----------------------------


sbatch timepoint_pcoa.sbatch ./csv_files/skin_floor_timepoint/lognorm_data_all.csv ./output/skin_floor_timepoint/PCOA/skin_floor \
"hand" floor 

#Rscript ./scripts/timepoint_pcoa.R ./csv_files/skin_floor_timepoint/lognorm_data_all.csv ./output/skin_floor_timepoint/PCOA/skin_floor \
#"hand" floor &> "./logs/PERMANOVA_PCOA/PERMANOVA_PCOA_timepoint_hand_floor.log"


## -------------------------------------

sbatch timepoint_pcoa.sbatch ./csv_files/skinassoc_timepoint/lognorm_data_all.csv ./output/skinassoc_timepoint/PCOA/skinassoc \
hand "hand associated"

#Rscript ./scripts/timepoint_pcoa.R ./csv_files/skinassoc_timepoint/lognorm_data_all.csv ./output/skinassoc_timepoint/PCOA/skinassoc \
#hand "hand associated" &> "./logs/PERMANOVA_PCOA/PERMANOVA_PCOA_timepoint_hand_hand associated.log"