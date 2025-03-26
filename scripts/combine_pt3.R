suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
source("./scripts/combine_pt3_utils.R")

args <- commandArgs(trailingOnly = TRUE)
folder = args[1]
pheno1 = args[2]
pheno2 = args[3]
NA_vals = args[4]
`%ni%` <- Negate(`%in%`)

#### open reference tables #####
combine_otus <- read.csv("./csv_files/combine/combine_otus.csv", check.names=FALSE) 
ontology <- read.csv("./csv_files/common_ontology.csv")

# assign NA as 0 if NA_vals is NOT provided 
if(is.na(NA_vals)){
combine_otus[is.na(combine_otus)] <- 0
}

ref.files <- list.files(".", "*_SraRunTable.txt", recursive = TRUE)
print("metadata files:")
print(ref.files)
print("_________________________")

# # PRJEB14474 (10172): 
# details_PRJEB14474 <- info(1, ref.files)
# combine_otus <- columns_to_filter(combine_otus, c("sample_name", "Date"), details_PRJEB14474, "2016-08-15 09:00")


# file PRJEB26708 (11470):
details_PRJEB26708 <- info(2,ref.files)
# combine_otus <- columns_to_filter(combine_otus, c("sample_name", "collection_timestamp"), details_PRJEB26708, "2016-09-11 09:00")
combine_otus <- columns_to_filter(combine_otus, c("sample_name", "collection_timestamp"), details_PRJEB26708, "2016-08-15 09:00")

# # file PRJEB33050 (12470):
details_PRJEB33050 <- info(3, ref.files)
combine_otus <- columns_to_filter(combine_otus, c("sample_name", "timepoint"), details_PRJEB33050, 2)


# file PRJEB6292 (2192):
details_PRJEB6292 <- info(4, ref.files)
combine_otus <- columns_to_filter(combine_otus, c("sample_name", "Day"), details_PRJEB6292, "D02")




#### Step 1.5: (in combine_otus) switch rows and columns ####

rownames(combine_otus) <- combine_otus$Genus
combine_otus <- combine_otus %>% subset(select = -c(Genus)) %>% t() %>% as.data.frame()
combine_otus <- tibble::rownames_to_column(combine_otus, "sample_name")
combine_otus$Phenotype <- NA
combine_otus$Study_ID <- NA


print("step 1 done")
print("________________________________________________________")


#### files ####

### PRJEB14474 (10172): ###

details_PRJEB14474 <- info(1, ref.files)
# make count df
print("count of each sample type in PRJEB14474: ")
count(details_PRJEB14474, sample_type) %>% print()


details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, "PRJEB14474", details_PRJEB14474, "sample_type", c(pheno1, pheno2))


# k in phenotypes df must have the same name as the ontology df
phenotypes <- data.frame(k = c(pheno1, pheno2),
                         Phenotype = c(1, 0))

combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("sample_type", "sample_name"), 
	"PRJEB14474", phenotypes, details_df_and_unneeded_phenotypes$vector)


print("________________________________________________________")
### PRJEB26708 (11740): ###

## details_PRJEB26708 is made above

print("count of each sample type in PRJEB26708 (11740) at interval 2016-08-15 09:00: ")
filter(details_PRJEB26708, collection_timestamp == "2016-08-15 09:00") %>% count(surface_sampledsample) %>% print()


details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, "PRJEB26708", details_PRJEB26708 , "surface_sampledsample", c(pheno1, pheno2))

# k in phenotypes df must have the same name as the ontology df
phenotypes  <- data.frame(k = c(pheno1, pheno2),
                          Phenotype = c(1, 0))


combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("surface_sampledsample", "sample_name"), 
	"PRJEB26708", phenotypes, details_df_and_unneeded_phenotypes$vector)



print("________________________________________________________")
### PRJEB33050 (12470): ###
# details_PRJEB33050 <- info(3, ref.files)

print("count of each sample type in PRJEB33050 (12470): ")
count(details_PRJEB33050, surfaceaggregated) %>% print()


details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, "PRJEB33050", details_PRJEB33050, "surfaceaggregated", c(pheno1, pheno2))

# k in phenotypes df must have the same name as the ontology df
phenotypes  <- data.frame(k = c(pheno1, pheno2),
                          Phenotype = c(1, 0))

combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("surfaceaggregated", "sample_name"), 
	"PRJEB33050", phenotypes, details_df_and_unneeded_phenotypes$vector)



print("________________________________________________________")

### PRJEB6292 (2192): ###
## details_PRJEB6292 is made above
print("count of each sample type in PRJEB6292 (2192) at interval D02: ")
filter(details_PRJEB6292, Day == "D02") %>% count(surface) %>% print()

details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, "PRJEB6292", details_PRJEB6292, "surface", c(pheno1, pheno2))

# k in phenotypes df must have the same name as the ontology df
phenotypes  <- data.frame(k = c(pheno1, pheno2),
                          Phenotype = c(1, 0))
# print(details_df_and_unneeded_phenotypes$df)

combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("surface", "sample_name"), 
	"PRJEB6292", phenotypes, details_df_and_unneeded_phenotypes$vector)

print("________________________________________________________")


#### Fixing an issue specific to this dataset

print("Fixing an issue specific to this dataset")

meta <- read.csv("./PRJEB14474/PRJEB14474_SraRunTable_old.txt", sep=",")
combine_otus_no_na <- filter(combine_otus, is.na(Study_ID) | is.na(Phenotype))

print("NA samples not from PRJEB14474")
filter(combine_otus_no_na, !(sample_name %in% meta$Run)) %>% dim() %>% print()
print("so all the NA samples are from PRJEB14474. This is because I filtered PRJEB14474's SRA table prior to analysis")

PRJEB14474 <- filter(combine_otus, Study_ID == "PRJEB14474" & (sample_name %in% details_PRJEB14474$sample_name))
combine_otus <- filter(combine_otus, Study_ID != "PRJEB14474")
combine_otus <- rbind(PRJEB14474, combine_otus)

print("________________________________________________________")

#### move Study_ID and Phenotype ####

combine_otus <- combine_otus %>% relocate(Phenotype, .after = sample_name) %>% relocate(Study_ID, .after = Phenotype)

print("before filtering")
print("Study IDs")
print(distinct(combine_otus, Study_ID))
print("Phenotypes")
print(distinct(combine_otus, Phenotype))
combine_otus %>% group_by(Study_ID) %>% distinct(Phenotype) %>% print()
count(combine_otus, Study_ID)
combine_otus %>% group_by(Phenotype) %>% count(Study_ID)


print("________________________________________________________")
#### Remove bacteria with no observations ####
# now that filtering of samples occurred, remove columns with no samples:

final_df <- combine_otus

# remove any cols with all NA... this can happen to low count rows that happened to have all assigned Genus's assigned as NA by dblur
cols = colSums(final_df[,-c(1:3)], na.rm = T) > 0
meta = c(sample_name = T,Phenotype = T, Study_ID = T)
final_df = final_df[c(meta, cols)]


# remove any rows with all NA... this can happen to low count rows that happened to have all assigned Genus's assigned as NA by dblur
# final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(x==0)),]
# final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(is.na(x))),]

rows = rowSums(final_df[,-c(1:3)], na.rm = T) > 0
final_df = final_df[rows,]


print("after filtering")
print("Study IDs")
print(distinct(final_df, Study_ID))
print("Phenotypes")
print(distinct(final_df, Phenotype))
final_df %>% group_by(Study_ID) %>% distinct(Phenotype) %>% print()
count(final_df, Study_ID)
final_df %>% group_by(Phenotype) %>% count(Study_ID)


print("PRJEB14474")
sample_count(final_df, details_PRJEB14474, "sample_type")

## details_PRJEB26708 is filter to the correct timepoint by details_df_and_unneeded_phenotypes
print("               ")
print("PRJEB26708")
sample_count(final_df, details_PRJEB26708, "surface_sampledsample")

print("               ")
print("PRJEB33050")
sample_count(final_df, details_PRJEB33050, "surfaceaggregated")

## details_PRJEB6292 is filter to the correct timepoint by details_df_and_unneeded_phenotypes
print("               ")
print("PRJEB6292")
sample_count(final_df, details_PRJEB6292, "surface")



# write.csv(final_df, "./csv_files/combine/final_df.csv", row.names = FALSE)



#### write output to file ####



lognorm(final_df[4:length(final_df)], final_df, paste("./csv_files/",folder,"/lognorm_data.csv",sep=""))

print("script complete")
