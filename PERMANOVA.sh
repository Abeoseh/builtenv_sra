#!/bin/bash

## just so I don't have to keep rewriting the commands everytime

sbatch PERMANOVA_PCOA.sbatch ./csv_files/associated_longitudinal/lognorm_data.csv \
./output/associated_longitudinal/PERMANOVA_PCOA/associated \
"hand associated" "floor"

sbatch PERMANOVA_PCOA.sbatch ./csv_files/skin_floor_longitudinal/lognorm_data.csv \
./output/skin_floor_longitudinal/PERMANOVA_PCOA/skin_floor \
"hand" "floor"

sbatch PERMANOVA_PCOA.sbatch ./csv_files/skinassoc_longitudinal/lognorm_data.csv \
./output/skinassoc_longitudinal/PERMANOVA_PCOA/skin_skinassoc \
"hand" "hand associated"


sbatch PERMANOVA.sbatch ./csv_files/pcoa_longitudinal/lognorm_data.csv \
./output/PCOA_longitudinal/PCOA_No_PRJEB14474 \
"hand" "environmental"







sbatch timepoint_pcoa.sbatch ./csv_files/associated_timepoint/lognorm_data_all.csv \
./output/associated_timepoint/PCOA/associated "hand associated" floor


sbatch timepoint_pcoa.sbatch ./csv_files/skin_floor_timepoint/lognorm_data_all.csv \
./output/skin_floor_timepoint/PCOA/skin_floor "hand" floor

sbatch timepoint_pcoa.sbatch ./csv_files/skinassoc_timepoint/lognorm_data_all.csv \
./output/skinassoc_timepoint/PCOA/skinassoc hand "hand associated"


sbatch timepoint_pcoa.sbatch ./csv_files/skinassoc_timepoint/lognorm_data.csv \
./output/skinassoc_timepoint/PCOA_RandomForestDataset/skinassoc hand "hand associated" T

sbatch timepoint_pcoa.sbatch ./csv_files/associated_timepoint/lognorm_data.csv \
./output/associated_timepoint/PCOA_RandomForestDataset/associated "hand associated" floor

sbatch timepoint_pcoa.sbatch ./csv_files/skin_floor_timepoint/lognorm_data.csv \
./output/skin_floor_timepoint/PCOA_RandomForestDataset/skin_floor "hand" floor



#### visualizing timepoints that weren't predicted well well after random forest ####
sbatch timepoint_pcoa.sbatch ./csv_files/associated_timepoint/lognorm_data.csv \
./output/associated_timepoint/PCOA_RandomForestDataset/associated "hand associated" floor T


sbatch timepoint_pcoa.sbatch ./csv_files/skin_floor_timepoint/lognorm_data.csv \
./output/skin_floor_timepoint/PCOA_RandomForestDataset/skin_floor "hand" floor T

sbatch timepoint_pcoa.sbatch ./csv_files/skinassoc_timepoint/lognorm_data.csv \
./output/skinassoc_timepoint/PCOA_RandomForestDataset/skinassoc hand "hand associated" T







