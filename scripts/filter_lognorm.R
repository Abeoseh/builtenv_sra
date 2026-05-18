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

input = "./csv_files/skin_floor_timepoint/lognorm_data_all.csv"
output = "./csv_files/skin_floor_timepoint/lognorm_data.csv"
# for_removal = c("PRJEB26708", "PRJEB6292") ## if you uncomment this, uncomment line 23

# print("removed datasets:")
# print(for_removal)

lognorm <- read.csv(input, check.names=FALSE)

# combine_otus <- filter(lognorm, !(Study_ID %in% for_removal)) ## if you uncomment this, uncomment line 16

## alternate use for this script, to filter time points without at least 30 samples
phenotype_counts = group_by(lognorm, Phenotype, timepoint, Study_ID) %>% count(Phenotype)
write.csv(phenotype_counts, "phenotype_counts.csv", row.names = F)
timepoints_for_removal = phenotype_counts$timepoint[phenotype_counts$n < 30]
print("timepoints that will be removed:")
print(timepoints_for_removal)
combine_otus <- lognorm[!lognorm$timepoint %in% timepoints_for_removal,]
print(colnames(combine_otus)[1:7])


print("before filtering")
print("Study IDs")
print(distinct(combine_otus, Study_ID))
print("Phenotypes")
print(distinct(combine_otus, Phenotype))
combine_otus %>% group_by(Study_ID, timepoint) %>% distinct(Phenotype) %>% print()
combine_otus %>% group_by(Study_ID) %>% distinct(Phenotype) %>% print()

final_df <- combine_otus

if ("timepoint" %in% colnames(final_df)){
  input = c(unlist(strsplit(input, "/")))[3]
  final_df %>% group_by(timepoint, Study_ID) %>% count(Phenotype) %>% write.csv(paste("output/",input,"/counts_after_filtering.csv",sep=""), row.names=F)
  print("entered and filtering")
  cols = colSums(final_df[,-c(1:6)], na.rm = T) > 0

  
  meta = c(sample_name = T,Phenotype = T, Study_ID = T, case = T, ID = T, timepoint = T)
  final_df = final_df[c(meta, cols)]
  print(colnames(final_df)[1:8])
  
  # remove any rows with all NA... this can happen to low count rows that happened to have all assigned Genus's assigned as NA by dblur
  # final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(x==0)),]
  # final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(is.na(x))),]
  
  rows = rowSums(final_df[,-c(1:6)], na.rm = T) > 0
  final_df = final_df[rows,]

}else{
  cols = colSums(final_df[,-c(1:5)], na.rm = T) > 0
  meta = c(sample_name = T,Phenotype = T, Study_ID = T)
  final_df = final_df[c(meta, cols)]
  
  
  # remove any rows with all NA... this can happen to low count rows that happened to have all assigned Genus's assigned as NA by dblur
  # final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(x==0)),]
  # final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(is.na(x))),]
  
  rows = rowSums(final_df[,-c(1:5)], na.rm = T) > 0
  final_df = final_df[rows,]
}


# final_df <- combine_otus
# filtered_cols = c()
# filtered_indicies = c()
# 
# for(i in 7:ncol(combine_otus)){
#   if(sum(as.array(combine_otus[[colnames(combine_otus)[i]]]), na.rm = TRUE) <= 0){
#     filtered_cols <- append(filtered_cols, colnames(combine_otus)[i])
#     filtered_indicies <- append(filtered_indicies, i)
#   }
# 
# 
# }
# 
# 
# filtered_cols <- data.frame("columns" = filtered_cols, "indices" = filtered_indicies)
# 
# final_df <- final_df[,-c(filtered_indicies)]
# 
# 
# # remove any rows with all NA... this can happen to low count rows that happened to have all assigned Genus's assigned as NA by dblur
# final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(x==0)),]
# final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(is.na(x))),]

print("______________________________________________________________________________")

print("after filtering")
print("Study IDs")
print(distinct(final_df, Study_ID))
print("Phenotypes")
print(distinct(final_df, Phenotype))
final_df %>% group_by(Study_ID) %>% distinct(Phenotype) %>% print()


count(final_df, Study_ID)

final_df %>% group_by(Phenotype) %>% count(Study_ID)

write.csv(final_df, "test.csv", row.names = F)

IDs <- unique(final_df$Study_ID)
if (length(IDs) != length(unique(lognorm$Study_ID))){
  print("adjusting IDs")
  final_df$ID = 0
  for(i in 1:length(IDs)){
      
  
      # print(distinct(final_df, ID))
      print(i)
      final_df$ID[final_df$Study_ID == IDs[i]] <- i-1
  }
}



write.csv(final_df, output, row.names=F)
print("script complete")

