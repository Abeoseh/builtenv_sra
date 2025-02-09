#!/usr/bin/env/Rscript

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))


### this is what happens when I do bind_rows(combine_otus, df) ###

# s a b c 
# x 1 2 3

# s d a c  e 
# y 4 5 2  6
# x 2 1 NA NA

# s a b  c  d  e
# x 1 2  3  NA NA
# y 5 NA 2 4  6
# x 1 NA NA 2  NA

files = list.files("./csv_files/combine", "single_", full.name=TRUE)
print("all files")
print(files)


combine_otus = read.csv(files[1], check.names=FALSE)

for(file in files[2:length(files)]){
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("current file")
	print(file)
	current_file = read.csv(file, check.names=FALSE)
	# check if all but the first column is in combine_otus 
	duplicate_samples <- (names(current_file)[-1] %in% names(combine_otus))
	# get the duplicate genera
	duplicate_genus <- current_file$Genus %in% combine_otus$Genus
	
	if( (TRUE %in% duplicate_samples )){
		print("duplicate samples exist")
		# index the names of the duplicates
		print(names(current_file)[-1][duplicate_samples])

		# add false for Genus column
		duplicate_samples = c(FALSE, duplicate_samples)
		# subset the current file by duplicate_samples so only NON-duplicate samples remain using !
		current_file <- current_file[,!duplicate_samples]
		combine_otus <- merge(combine_otus, current_file, by="Genus", all=TRUE)
		
		# old method
		#combine_otus <- bind_rows(combine_otus, current_file) 
		#if( length(unique(combine_otus$Genus)) != nrow(combine_otus) ){
			#print("duplicate genus exist!!")
			# to print the duplicate Genus
			# n_occur <- data.frame(table(combine_otus$Genus))
			# print( combine_otus$Genus[combine_otus$Genus %in% n_occur$Var1[n_occur$Freq > 1]] )
			# write.csv(combine_otus, "./logs/combine/error.csv", row.names=FALSE)
			#combine_otus <- aggregate(. ~ Genus, combine_otus, sum, na.rm=TRUE, na.action=NULL)}
	}

	else( combine_otus <- merge(combine_otus, current_file, by="Genus", all=TRUE) )
}

# print("combine df")
# print(colnames(combine_otus))
# print(head(combine_otus, 10))
print("number of genera and unique genera")
print(length(unique(combine_otus$Genus)))
print(nrow(combine_otus))

write.csv(combine_otus, "./csv_files/combine/combine_otus.csv", row.names=FALSE)
print("done")
