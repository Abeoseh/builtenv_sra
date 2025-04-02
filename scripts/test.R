suppressPackageStartupMessages(library(dplyr))
setwd("C:/Users/brean/Downloads/masters/Fodor/builtenv_sra/")

df <- read.csv("./csv_files/comparisons.csv")


intersect(df$feature_dorm, df$feature_hospital)


df <- read.csv("./PRJEB26708/PRJEB26708_SraRunTable.csv")

# surface_sampledsample
df %>% count(collection_timestamp)
# View

df %>% group_by(collection_timestamp) %>% count(surface_sampledsample) %>% View()
# filter(collection_timestamp == "2016-09-11 09:00" | collection_timestamp == "2016-08-15 09:00") %>% 

?intersect


# PRJEB14474
# PRJEB6292
auc.df <- read.csv("./output/skin_floor/CSVs/builtenv_AUCs.csv") %>% filter(Study_ID == "PRJEB6292")

a <- auc.df[auc.df$Permutation == FALSE,]$AUC
samp <- auc.df[auc.df$Permutation == TRUE,]$AUC
length(samp[samp >= a])/1000  

length(samp[samp > a]) + length(samp[samp < a]) + length(samp[samp == a])
