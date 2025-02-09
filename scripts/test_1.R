#!/usr/bin/env/Rscript

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(biomformat))
# suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))


# DEBIAS_data <- read.csv("./output/sink_nonsink/DEBIAS-M_runs/builtenv_debiased_lognorm_PRJNA878661.csv", colClasses = c("Phenotype" = "factor"))
# 
# DEBIAS_data %>% distinct(Phenotype) %>% print()


new <- read.csv("./PRJEB14474/PRJEB14474_SraRunTable.txt", sep=",")
# old <- read.csv("./csv_files/skin_floor/lognorm_data.csv",sep=",")

new$Date <- as.Date(new$Date, format = "%M/%D/%Y")

new <- filter(new, Date > 2/23/2013)
write.csv(new, "csv_files/test/test.csv", row.names=F)


# figure out which samples are missing
# df <- filter(df, Day == "D02") # %>% filter(surface %in% c("Bedroom_Floor",	"Kitchen_Floor",	"Bathroom_Floor"))

# IDs <- df$Run[ !(df$Run %in% df1$sample_name) ]
# 
# df %>% filter(Run %in% IDs) %>% select(Run, surface)



# df <- read.csv("./csv_files/sink_nonsink/lognorm_data.csv")
# dim(df) %>% print()
# count(df, Study_ID)
# 
# df %>% group_by(Phenotype) %>% count(Study_ID)

# l <- rowSums(df[,6:length(df)])
# print(l[l == 0])

# rows_with_sum_zero <- rowSums(df[,6:length(df)]) == 0
# df[rows_with_sum_zero, ] %>% print()