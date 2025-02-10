#!/usr/bin/env/Rscript

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(biomformat))
#suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]

print(paste("input folder:",input))
print(paste("output folder:",output))


lognorm <- read.csv(paste("./csv_files/",input,"/lognorm_data.csv",sep=""), check.names=FALSE)
phenos <- read.csv("./csv_files/phenotypes.csv") # for naming the graphs
IDs = distinct(lognorm, Study_ID)$Study_ID


## Aims:
# 1: preform all the wilcox between skin and floor counts for one study
# 2: graph it against the p-value for the wilcox of the other studies

compute_pvals <- function(subsetted_data, Study_ID){
	
	# create wilcox_pval df
	wilcox_pvals = data.frame(matrix(nrow = 0, ncol = 4))
	colnames(wilcox_pvals) <- c("bacteria", paste("pval_",Study_ID, sep=""), paste("skin_mean_",Study_ID, sep=""), paste("floor_mean_",Study_ID, sep=""))
	
	# bacteria that exists within the Study_ID
	all_cols = colnames(subsetted_data)
	bacteria_cols = all_cols[6:length(all_cols)]

	# skin is 1, floor is 0
	skin <- filter(subsetted_data, Phenotype == 1)
	floor <- filter(subsetted_data, Phenotype == 0)
	

	for( bacteria in bacteria_cols ){
		# the !! unpacks the variable bacteria so R is not literally using "bacteria"
		# make a vector of sample counts for skin and floor
		current_skin <- select(skin, !!bacteria)[[bacteria]]
		current_floor  <- select(floor, !!bacteria)[[bacteria]]
		#print("dimensions after selecting")
		#print(sum(current_skin))
		#print(sum(current_floor))
	
		# compute the mean for skin and floor without NA values
		floor_mean = mean(current_floor, na.rm = T)
		skin_mean = mean(current_skin, na.rm = T)
		# conditionally change p-values to be positive or negative depending on the means
		# compute mean of current_skin and current_floor, if mean in current_skin is higher log10(pval) else, -log10(pval)
		if( is.na(floor_mean) | is.na(skin_mean) ){
			#print("ID bacteria floor_mean skin_mean")
			#print(paste(Study_ID, bacteria, floor_mean, skin_mean))

			wilcox_pvals[nrow(wilcox_pvals) +1,] <- c( bacteria, NaN, NaN, NaN ) }
		else{
			pval = wilcox.test(current_skin, current_floor, alternative = "two.sided", paired = FALSE)$p.value 

			if (is.na(pval)){
				wilcox_pvals[nrow(wilcox_pvals) +1,] <- c( bacteria, NaN, skin_mean, floor_mean ) }

			else if( skin_mean > floor_mean ){
				wilcox_pvals[nrow(wilcox_pvals) +1,] <- c( bacteria, -log10(pval), skin_mean, floor_mean ) }

			else{ wilcox_pvals[nrow(wilcox_pvals) +1,] <- c( bacteria, log10(pval), skin_mean, floor_mean ) }
		}
	}
	return(wilcox_pvals)
}



print("IDs")
print(IDs)

### 10172 ### 
# https://stackoverflow.com/questions/34219912/how-to-use-a-variable-in-dplyrfilter
# filter was trying to literally compare ID[1] which doesn't exist... instead I needed to "unpack" ID[1] to 10172 so filter used that instead

wilcox_pval =  filter(lognorm, Study_ID == !!IDs[1]) %>% compute_pvals(IDs[1])

print(paste("done", IDs[1]))

for( ID in IDs[2:length(IDs)]){
	wilcox_pval = filter(lognorm, Study_ID == !!ID) %>% compute_pvals(ID) %>% merge(wilcox_pval, by="bacteria", all=TRUE)
	print(paste("done", ID))
}

# warnings()

## clean up
print("str of wilcox_pval before numeric conversion")
print(str(wilcox_pval))
wilcox_pval[,2:(length(wilcox_pval))] <- sapply(wilcox_pval[,2:(length(wilcox_pval))], as.numeric)
print("str of wilcox_pval after numeric conversion (NAs are expected)")
print(str(wilcox_pval))

# write.csv(wilcox_pval, "./csv_files/test/wilcox_pval.csv", row.names = FALSE)
write.csv(wilcox_pval, paste("./output/", output, "/pval_v_pval/files/wilcox_pval.csv", sep=""),row.names=F )


### Graph: ###
used_IDs = list()
for( ID1 in IDs ){
	used_IDs = append(used_IDs, ID1)

	for( ID2 in IDs ){
		if( ID1 != ID2 & !(ID2 %in% used_IDs)){
			print(paste(ID1, "vs", ID2))
			pval_ID1 = paste("pval_",ID1,sep="")
			pval_ID2 = paste("pval_",ID2,sep="")
			df <- select(wilcox_pval, bacteria, !!pval_ID1, !!pval_ID2) 
			df <- na.omit(df)

			# make the x and y ranges the same
			mylims <- range(c(df[[pval_ID1]], df[[pval_ID2]]))
			print("x and y limits")
			print(mylims)
			
			# for naming the axes 
			phen1 <- filter(phenos, ID == ID1)
			phen1 <- phen1[1,2]

			phen2 <- filter(phenos, ID == ID2)
			phen2 <- phen2[1,2]


			png(paste("./output/",output,"/pval_v_pval/plot_",ID1,"v",ID2,".png",sep=""), width = 1050, height = 480)

			plot = ggplot(df, aes(x = df[[pval_ID1]], y = df[[pval_ID2]], label = bacteria)) +
				geom_point() +
				# geom_text(mapping = aes(label = bacteria)) + # remove label from ggplot() aes before you uncomment this line of code (2 lines up)
				geom_hline(yintercept=log10(0.05), linetype='dotted', col = 'red') +
				geom_hline(yintercept=-log10(0.05), linetype='dotted', col = 'red') +
				geom_vline(xintercept=log10(0.05), linetype='dotted', col = 'red') +
				geom_vline(xintercept=-log10(0.05), linetype='dotted', col = 'red') +
				# make the x and y ranges the same
				# coord_cartesian(xlim = mylims, ylim = mylims) +
				geom_text_repel(max.overlaps = 10, force_pull = 1, nudge_y = 1,size = 3) +
				labs(title = "log10 p-value vs p-value plot", x = phen1, y = phen2) +
				theme(plot.title = element_text(size=22), axis.text=element_text(size=11),
        			axis.title=element_text(size=15)) 
				

			print(plot)

			dev.off()
		}
	}
}

print("done")


