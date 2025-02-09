#!/usr/bin/env/Rscript

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(biomformat))
# suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(vegan))

print("start")

input = "./csv_files/associated/lognorm_data_all.csv"
output = "./csv_files/associated/lognorm_data.csv"
for_removal = c("PRJEB26708", "PRJEB6292")

print("removed datasets:")
print(for_removal)

lognorm <- read.csv(input, check.names=FALSE)

combine_otus <- filter(lognorm, !(Study_ID %in% for_removal))


print("before filtering")
print("Study IDs")
print(distinct(lognorm, Study_ID))
print("Phenotypes")
print(distinct(lognorm, Phenotype))
lognorm %>% group_by(Study_ID) %>% distinct(Phenotype) %>% print()



final_df <- combine_otus
filtered_cols = c()
filtered_indicies = c()

for(i in 4:ncol(combine_otus)){
  if(sum(as.array(combine_otus[[colnames(combine_otus)[i]]]), na.rm = TRUE) <= 0){
    filtered_cols <- append(filtered_cols, colnames(combine_otus)[i])
    filtered_indicies <- append(filtered_indicies, i)
  }


}


filtered_cols <- data.frame("columns" = filtered_cols, "indices" = filtered_indicies)

final_df <- final_df[,-c(filtered_indicies)]


# remove any rows with all NA... this can happen to low count rows that happened to have all assigned Genus's assigned as NA by dblur
final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(x==0)),]
final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(is.na(x))),]


print("after filtering")
print("Study IDs")
print(distinct(final_df, Study_ID))
print("Phenotypes")
print(distinct(final_df, Phenotype))
final_df %>% group_by(Study_ID) %>% distinct(Phenotype) %>% print()


count(final_df, Study_ID)

final_df %>% group_by(Phenotype) %>% count(Study_ID)

IDs <- unique(final_df$Study_ID)
final_df$ID = 0
for(i in 1:length(IDs)){

    # print(distinct(final_df, ID))
    print(i)
    final_df$ID[final_df$Study_ID == IDs[i]] <- i-1
}




write.csv(final_df, output, row.names=F)

# DEBIAS_data <- read.csv("./output/sink_nonsink/DEBIAS-M_runs/builtenv_debiased_lognorm_PRJNA878661.csv", colClasses = c("Phenotype" = "factor"))
#
# DEBIAS_data %>% distinct(Phenotype) %>% print()


## determining which samples from the ontology didn't make it 

# df <- read.csv("./PRJEB14474/PRJEB14474_SraRunTable.txt", sep=",")
# df1 <- read.csv("./csv_files/sink_nonsink/lognorm_data.csv")
# 
# df <- filter(df, sample_type %in% c("Cold Tap Faucet Handle", "Countertop", "Hot Tap Water", "Hot Tap Water Faucet Handle",	"Cold Tap Water Faucet Handle",	"Cold Tap Water"))
# dim(df)
# df1 <- filter(df1, Study_ID == "PRJEB14474")
# dim(df1)
# 
# 
# IDs <- df$Run[ !(df$Run %in% df1$sample_name) ]
# 
# 
# 
# df %>% filter(Run %in% IDs) %>% select(Run, sample_type)
# 
# 
## study breakdown by group
# df <- read.csv("./csv_files/sink_nonsink/lognorm_data.csv")
#
# dim(df) %>% print()
# count(df, Study_ID)
#
# df %>% group_by(Phenotype) %>% count(Study_ID)





# l <- rowSums(df[,6:length(df)])
# print(l[l == 0])

# rows_with_sum_zero <- rowSums(df[,6:length(df)]) == 0
# df[rows_with_sum_zero, ] %>% print()