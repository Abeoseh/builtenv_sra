suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
folder = args[1]
pheno1 = args[2]
pheno2 = args[3]
NA_vals = args[4]
`%ni%` <- Negate(`%in%`)

#### open reference tables #####

combine_otus <- read.csv("./csv_files/combine/combine_otus.csv", check.names=FALSE) 

ontology <- read.csv("./csv_files/combine/common_ontology_2192.csv")

# assign NA as 0 if NA_vals is NOT provided 
if(is.na(NA_vals)){
combine_otus[is.na(combine_otus)] <- 0
}
ref.files <- list.files(pattern =".txt$", recursive = TRUE)



#### function 1-2: read the metadata files and removed undesired samples (rows) ####
info <- function(ref.file_num, display_cols=NULL){
  ### read the metadata files ###
  metadata <- read.csv(ref.files[ref.file_num], sep = ",", header=TRUE)
  colnames(metadata)[1] = "sample_name"
  
  if(!is.null(display_cols)){
    print("metadata columns:")
    print(colnames(metadata))
  }
  
  print(paste("ref_file:", ref.files[ref.file_num]))  
  
  return(metadata)
}



columns_to_filter <- function(df, columns, ref, condition){
### filter values out of columns based on provided condition ###
#  ref <- ref

  ref <- ref[,columns]

  r2 <- ref[ref[[columns[2]]] != condition,]
    
  r <- df[,!(names(df) %in% r2[[columns[1]]])]

  return(r)

}

# file 11470:

details_11740 <- info(2)
combine_otus <- columns_to_filter(combine_otus, c("sample_name", "collection_timestamp"), details_11740, "2016-08-15 09:00")
print("done all")

# file 2192
details_2192<- info(4)
combine_otus <- columns_to_filter(combine_otus, c("sample_name", "day"), details_2192, "D02")

#### Step 1.5: (in combine_otus) switch rows and columns ####

rownames(combine_otus) <- combine_otus$Genus
combine_otus <- combine_otus %>% subset(select = -c(Genus)) %>% t() %>% as.data.frame()
combine_otus <- tibble::rownames_to_column(combine_otus, "sample_name")
combine_otus$Phenotype <- NA
combine_otus$Study_ID <- NA


print("step 1 done")
#print("head of combine_otus")
#print(head(combine_otus))
print("________________________________________________________")

##### function 3-4: Add a common ontology and merge to combine df ####

unneeded_phenotypes <- function(ontology_df, current_ID, details_df, details_phenotype, ontology){

	# ontology_df is a dataframe with the following columns: 
	### ID: contains all the Study_IDs
	### surface: the original sample names ex bed rail, kitchen table, floor corner
	### common_name: the ontology between studies (a.k.a chosen common names)... all phenotypes that are unneeded are NA!! 
	
	# current_ID is the Study ID of the current study
	# details_df is the metadata file
	# details_phenotype is the name of the phenotype column in details_df
	# ontology is a vector with the desired ontology

	# example of how to run
	### details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, 10172, details_10172, "sample_type", c("skin associated", "floor associated"))
    print("ontology")
    print(ontology)

	current_ontology <- select(ontology_df, ID, surface, common_name) %>% filter(ID == current_ID & !(is.na(common_name))) 

	# select the rows of current_ontology which have values labeled as your first ontology phenotype and...
	pheno_1_associated <- current_ontology$surface[current_ontology$common_name == ontology[1]]

	# rename those columns in details_df as that ontology.
	details_df[[details_phenotype]][details_df[[details_phenotype]] %in% pheno_1_associated] <- ontology[1]
	
	# Do the same for the second phenotype.
	pheno_2_associated <- current_ontology$surface[current_ontology$common_name == ontology[2]]
	details_df[[details_phenotype]][details_df[[details_phenotype]] %in% pheno_2_associated] <- ontology[2]
	
	# Select the samples that are not part of the overall ontology
	uneeded_phenotypes = filter(details_df, .data[[details_phenotype]] %ni% ontology ) %>% select(!!details_phenotype) %>% unique()
	uneeded_phenotypes = uneeded_phenotypes[[details_phenotype]]
	print(paste("uneeded phenotypes",current_ID))
	print(uneeded_phenotypes)


	return(list(df = details_df, vector = uneeded_phenotypes))

}

add_info_cols <- function(r, ref, col_names, study_ID, ref.df, undesirable = NULL){

  # r: merged count files (combine_otus)
  # ref: individual metadata file
  # col_names: a vector of two column names from the metadata file containing: "sample_type", "sample_name"
  # study_ID: the study ID (either self assigned or from qitta)
  # ref.df: a data frame aligning the binary classes in sample_type with 0s and 1s. 
  # undesirable: a vector containing all undesirable classes from sample_type

  
  # read file
  ref <- ref
  print("done 1")
  # subset columns need: sample_Id, phenotype
  ref <- ref[,col_names]
  print("done 2")
  # only select the sample metadata I have data on. 
  
  ref <- ref[ref[[col_names[2]]] %in% r[[col_names[2]]], ]
  print("done 2")
  
  # remove undesirable samples
  if(!is.null(undesirable)){
    
    #temp <- ref[ref[[col_names[1]]] == undesirable,]
    
    #ref <- ref[ref[[col_names[1]]] != undesirable,]
    #r <- r[!r[["sample_name"]] %in% temp[[col_names[2]]], ]
    
    temp <- dplyr::filter(ref, .data[[col_names[1]]] %in% undesirable)
    print("done 3")
    ref <-  dplyr::filter(ref, .data[[col_names[1]]] %ni% undesirable)
    print("done 4")
    r <- r[!r[["sample_name"]] %in% temp[[col_names[2]]], ]
    print("done 5")
    print("ref values")
    print(distinct(ref, .data[[col_names[1]]]))
    
    
  }
  
  #add study ID column
  ref$Study_ID <- study_ID
  print("done 6")
  
  
  #  change reference column to 1s and 0s
  ref <- merge(ref, ref.df, by.x = col_names[1], by.y = "k", all.x=TRUE)
  print("done 7")
  
  
  ref <- ref[, 2:length(ref)]
  print("done 8")
  

  
  # merge ref to r
  r$Phenotype[match(ref$sample_name, r$sample_name)] <- ref$Phenotype
  r$Study_ID[match(ref$sample_name, r$sample_name)] <- ref$Study_ID

  ## ALT:
  #inters <- intersect(names(r), names(ref))
  #r <- merge(r, ref, by = inters, all = TRUE) 
  # r <- aggregate(.~sample_name, data= r, FUN=sum, na.rm=TRUE, na.action = NULL) 
  
  return(r)
}




#### files ####

### 10172:
details_10172 <- info(1)
print("count of each sample type in 10172: ")
count(details_10172, sample_type) %>% print()

# make count df
count = count(details_10172, sample_type)
count$Study_ID = 10172


details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, 10172, details_10172, "sample_type", c(pheno1, pheno2))


current_ontology <- select(ontology, ID, surface, common_name) %>% filter(ID == 10172 & !(is.na(common_name))) 


# k in phenotypes df must have the same name as the ontology df
phenotypes <- data.frame(k = c(pheno1, pheno2),
                         Phenotype = c(1, 0))

combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("sample_type", "sample_name"), 
	10172, phenotypes, details_df_and_unneeded_phenotypes$vector)

print("________________________________________________________")
### 11740: ###

print("count of each sample type in 11740: ")
count(details_11740, surface_sampledsample) %>% print()


current_ontology <- select(ontology, ID, surface, common_name) %>% filter(ID == 11740 & !(is.na(common_name))) 

details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, 11740, details_11740, "surface_sampledsample", c(pheno1, pheno2))

# k in phenotypes df must have the same name as the ontology df
phenotypes  <- data.frame(k = c(pheno1, pheno2),
                          Phenotype = c(1, 0))


combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("surface_sampledsample", "sample_name"), 
	11740, phenotypes, details_df_and_unneeded_phenotypes$vector)

print("________________________________________________________")
### 12470: ###
details_12470 <- info(3)

print("count of each sample type in 12470: ")
count(details_12470, description) %>% print()


current_ontology <- select(ontology, ID, surface, common_name) %>% filter(ID == 12470 & !(is.na(common_name))) 

details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, 12470, details_12470, "description", c(pheno1, pheno2))

# k in phenotypes df must have the same name as the ontology df
phenotypes  <- data.frame(k = c(pheno1, pheno2),
                          Phenotype = c(1, 0))

combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("description", "sample_name"), 
	12470, phenotypes, details_df_and_unneeded_phenotypes$vector)



print("________________________________________________________")

# 2192:
print("count of each sample type in 2192: ")
# details_2192 is made above
count(details_2192, subject_surface) %>% print()

current_ontology <- select(ontology, ID, surface, common_name) %>% filter(ID == 2192 & !(is.na(common_name))) 

details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, 2192, details_2192, "subject_surface", c(pheno1, pheno2))

# k in phenotypes df must have the same name as the ontology df
phenotypes  <- data.frame(k = c(pheno1, pheno2),
                          Phenotype = c(1, 0))

combine_otus <- add_info_cols(combine_otus, details_df_and_unneeded_phenotypes$df, c("subject_surface", "sample_name"), 
	2192, phenotypes, details_df_and_unneeded_phenotypes$vector)

print("________________________________________________________")



#write.csv(count, "./csv_files/combine/count.csv", row.names=FALSE)

#### move Study_ID and Phenotype ####

combine_otus <- combine_otus %>% relocate(Phenotype, .after = sample_name) %>% relocate(Study_ID, .after = Phenotype)

print("Study IDs")
print(distinct(combine_otus, Study_ID))
print("Phenotypes")
print(distinct(combine_otus, Phenotype))
combine_otus %>% group_by(Study_ID) %>% distinct(Phenotype) %>% print()


#write.csv(combine_otus, "./csv_files/combine/unfiltered_phenos_otus.csv", row.names = FALSE)




#### Remove bacteria with no observations ####
# now that filtering of samples occurred, remove columns with no samples:

print("________________________________________________________")

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


# remove any rows with all NA/0... this can happen to low count rows that happened to have all assigned Genus's assigned as NA by dblur
final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(x==0)),]
final_df <- final_df[apply(final_df[,-c(1:3)], 1, function(x) !all(is.na(x))),]


# write.csv(final_df, "./csv_files/combine/final_df.csv", row.names = FALSE)


#### write output to file ####


#### log normalizing and getting into DEBIAS-M format ####
lognorm <- function(table, dataframe, csv_file, filter = NULL, return_table = NULL)
{
  # actual lognorm
  avg <- sum(rowSums(table, na.rm = T))/nrow(table)
  table <- sweep(table,1,rowSums(table, na.rm = T),"/")
  table <- log10(table*avg + 1)

  #### filtering rare taxa #### all commented out since all the taxa kept being removed
  #print("dim lognorm before filtering rare taxa")
  #print(dim(table))
  #
  #percent = colSums(is.na(table))/nrow(table) # percent is a named number
  #percent = percent[percent < 0.25] # gives the columns with an amount of NAs less than 25%
  #
  #table <- table %>% select(names(percent)) # select the columns with the amount of NAs less than 25%
  #
  #print("dim lognorm after filtering NA")
  #print(dim(table))
  #
  #percent = colSums(table == 0, na.rm = T)#/nrow(table) # percent is a named number
  #percent = percent[percent < 0.25] # gives the columns with an amount of NAs less than 25%
  #write.csv(as.data.frame(percent), "./csv_files/test/percent.csv")
  #table <- table %>% select(names(percent)) # select the columns with the amount of 0s less than 25%

  #print("dim lognorm after filtering 0s")
  #print(dim(table))

  # add sample_name, Study_ID, Phenotype back   
  table <- add_column(table, Study_ID=dataframe$Study_ID, .before = colnames(table)[1])
  table <- add_column(table, Phenotype=dataframe$Phenotype, .before = colnames(table)[1])
  table <- add_column(table, sample_name=dataframe$sample_name, .before = colnames(table)[1])
  IDs <- distinct(table, Study_ID)$Study_ID

  # DEBIAS-M needs the study IDs to be labeled starting from 0 
  table$ID = 0
  for(i in 1:length(IDs)){
    
    # print(distinct(df, ID))
    # print(i)
    table$ID[table$Study_ID == IDs[i]] <- i-1
  }
  
  table$case <- case_when(
    table$Phenotype == 1 ~ TRUE,
    table$Phenotype == 0 ~ FALSE,
  )
  table <- relocate(table, ID, .after = Study_ID)
  table <- relocate(table, case, .after = Study_ID)
  
  # filter columns with all 0s
  if(!is.null(filter)){
    final_df <- table
    
    
    filtered_cols = c()
    filtered_indicies = c()

    for(i in 4:ncol(table)){
      if(sum(as.array(table[[colnames(table)[i]]]), na.rm = TRUE) <= 0){
        filtered_cols <- append(filtered_cols, colnames(table)[i])
        filtered_indicies <- append(filtered_indicies, i)
      } 
    }
    
    filtered_cols <- data.frame("columns" = filtered_cols, "indices" = filtered_indicies)
    final_df <- final_df[,-c(filtered_indicies)]
    write.csv(final_df, csv_file, row.names = FALSE)
  }
  
  else(write.csv(table, csv_file, row.names = FALSE))
  print(dim(table))
  if(!is.null(return_table)){return(table)}
  

}

lognorm(final_df[4:length(final_df)], final_df, paste("./csv_files/",folder,"/lognorm_data.csv",sep=""))

print("script complete")
