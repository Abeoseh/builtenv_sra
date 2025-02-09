#!/bin/bash

# just so I don't have to keep rewriting the commands everytime

sbatch PERMANOVA_PCOA.sbatch ./csv_files/associated/lognorm_data.csv \
./output/associated/PERMANOVA_PCOA/associated \
"skin associated" "floor associated"

sbatch PERMANOVA_PCOA.sbatch ./csv_files/skin_floor/lognorm_data.csv \
./output/skin_floor/PERMANOVA_PCOA/skin_floor \
"skin" "floor"

sbatch PERMANOVA_PCOA.sbatch ./csv_files/skin_skinassoc/lognorm_data.csv \
./output/skin_skinassoc/PERMANOVA_PCOA/skin_skinassoc \
"skin" "skin associated"
