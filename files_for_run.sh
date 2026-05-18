#!/bin/bash

## pval-pval plots
#sbatch pval-pval_plot.sbatch associated_na associated

#sbatch pval-pval_plot.sbatch skin_skinassoc_na skin_skinassoc


## random forest pre-DEBIAS
./run_all_pt1.sh skin_skinassoc skin_skinassoc 4

./run_all.sh associated associated 4

./run_all.sh skin_floor skin_floor 4

## PCoA and PERMANOVA
#sbatch PERMANOVA_PCOA.sbatch ./csv_files/skin_skinassoc/lognorm_data_all.csv \
#./output/skin_skinassoc/PERMANOVA_PCOA/skin_skinassoc \
#"hand" "hand associated"

#sbatch PERMANOVA_PCOA.sbatch ./csv_files/associated/lognorm_data_all.csv \
#./output/associated/PERMANOVA_PCOA/associated \
#"hand associated" "floor"

#sbatch PERMANOVA_PCOA.sbatch ./csv_files/skin_floor/lognorm_data.csv \
#./output/skin_floor/PERMANOVA_PCOA/skin_floor \
#"hand" "floor"


