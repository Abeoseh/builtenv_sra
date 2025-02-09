suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

# args <- commandArgs(trailingOnly = TRUE)
# folder = args[1]
# pheno1 = args[2]
# pheno2 = args[3]
# NA_vals = args[4]
`%ni%` <- Negate(`%in%`)

#### open reference tables #####

#### function 1-2: read the metadata files and removed undesired samples (rows) ####
info <- function(ref.file_num, ref.files, display_cols=NULL){
  ### read the metadata files ###
  metadata <- read.csv(ref.files[ref.file_num], sep = ",", header=TRUE) #, colClasses = c("sample_name" = "character"))
  colnames(metadata)[1] = "sample_name"
  
  if(!is.null(display_cols)){
    print("metadata columns:")
    print(colnames(metadata))
  }
  
  print(paste("ref_file:", ref.files[ref.file_num]))  
  
  return(metadata)
}


columns_to_filter <- function(df, columns, ref, condition){
### remove values from combine_otus columns based on provided condition ###
#  ref <- ref

  ref <- ref[,columns]

  r2 <- ref[ref[[columns[2]]] != condition,]
    
  r <- df[,!(names(df) %in% r2[[columns[1]]])]

  return(r)

}


##### functions 3&4: Add a common ontology and merge to combine df ####

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
	### details_df_and_unneeded_phenotypes <- unneeded_phenotypes(ontology, 10172, details_10172, "subject_surface", c("skin associated", "floor associated"))
  print("ontology")
  print(ontology)
  
	# current_ID <- deparse(substitute(current_ID))
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


sample_count <- function(data, metadata, pheno_col){
	## Small function which counts the amount of samples from the metadata file present in combine_outs
	# data is combine_otus
	# metadata is the SraRunTable
	# Study_ID is the study_ID within combine
	# pheno_col is the phenotype column within metadata
	
	# data <- filter(data, Study_ID == Study_ID) # maybe not needed
	# IDs <- metadata$sample_name[ (metadata$sample_name %in% data$sample_name) ] # maybe not needed
	metadata %>% filter(sample_name %in% data$sample_name) %>% select(sample_name, all_of(pheno_col)) %>% group_by(.data[[pheno_col]]) %>% count()
		


}


#### log normalizing and getting into DEBIAS-M format ####
lognorm <- function(table, dataframe, csv_file, filter = NULL, return_table = NULL){
  # actual lognorm
  avg <- sum(rowSums(table, na.rm = T))/nrow(table)
  table <- sweep(table,1,rowSums(table, na.rm = T),"/")
  table <- log10(table*avg + 1)

  # add sample_name, Study_ID, Phenotype back   
  table <- add_column(table, Study_ID=dataframe$Study_ID, .before = colnames(table)[1])
  table <- add_column(table, Phenotype=dataframe$Phenotype, .before = colnames(table)[1])
  table <- add_column(table, sample_name=dataframe$sample_name, .before = colnames(table)[1])
  IDs <- distinct(table, Study_ID)$Study_ID

  # DEBIAS-M needs the study IDs to be labeled starting from 0 
  table$ID = 0
  for(i in 1:length(IDs)){
    
    # print(distinct(df, ID))
    print(i)
    table$ID[table$Study_ID == IDs[i]] <- i-1
  }
  
  table$case <- case_when(
    table$Phenotype == 1 ~ TRUE,
    table$Phenotype == 0 ~ FALSE,
  )
  table <- relocate(table, ID, .after = Study_ID)
  table <- relocate(table, case, .after = Study_ID)
  
  write.csv(table, csv_file, row.names = FALSE)
  print("amount of rows and columns in df:")
  print(dim(table))
  if(!is.null(return_table)){return(table)}
  

}
