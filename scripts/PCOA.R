.libPaths( c( .libPaths(), "~/my_R_libs") )

library(dplyr)
library(stringr)
# PCoA
library(ecodist)
library(vegan)
# PCA
library(factoextra)
library(lemon)
library(ggplot2)

### IF statement explaination:
# the images are arranged in 2x2 grids

args = commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]

output = strsplit(output, ".png")[[1]]

PCOA_df <- read.csv(input)
Studies = read.csv("./csv_files/phenotypes.csv")
IDs = unique(PCOA_df$Study_ID)

pcoa <- function(df, chosen_title, pc = NULL){
  
  
  
  bray <- vegdist(df, method = "bray")
  pcoa_val <- pco(bray, negvals = "zero", dround = 0)
  
  if(!is.null(pc)){
    
    cat("example to access elements in the list: \n
    l = pcoa(df)\n
    l[1] # returns bray-curtis values\n
    l[2] # returns the pcoa results\n
    ls[2][[1]] # returns vectors and values which can be accessed via $\n
    ex: ls[2][[1]]$vectors")
    
    return(list(bray, pcoa_val))
  }
  
  pcoa_val.df = data.frame(sample_name = row.names(pcoa_val$vectors),
                           PCoA1 = pcoa_val$vectors[,1], 
                           PCoA2 = pcoa_val$vectors[,2])
  
  pcoa_val.df = merge(pcoa_val.df, phenos, by.all = "sample_name")
  
  # compute variance explained by each PCoA via eigenvalues 
  eigenvalues = pcoa_val$values
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = as.factor(Phenotype)), alpha = 0.7) +
    scale_color_brewer(name = "phenotype", palette = "Accent") +
	theme(plot.margin = margin(10, 10, 20, 10)) +
    labs(title = str_wrap(chosen_title,30), x = paste("PC1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
        y = paste("PC2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""))
  return(pco.plot)
}


for(i in 1:length(IDs)){
  post_DEBIAS <- filter(PCOA_df, Study_ID == IDs[i])

  
  phenos <- select(post_DEBIAS, sample_name, Phenotype)
  phenos$Phenotype[phenos$Phenotype=="1"] <- "sink"
  phenos$Phenotype[phenos$Phenotype=="0"] <- "non sink"
  
  row.names(post_DEBIAS) <- post_DEBIAS$sample_name
  post_DEBIAS <- post_DEBIAS[,6:length(post_DEBIAS)]
  row.names(post_DEBIAS) <- phenos$sample_name
  
  
  Study = Studies$Author[Studies$ID == IDs[i]]
  print(Study)

  # ex i==1, 1+1%%2 == 0 and 1+1%%4 is not 0
  if((i+1)%%2 == 0 & (i+1)%%4 != 0){
    p1 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
    print(paste(i, "p1"))
    amount = 1}
  
  if(i%%2 == 0 & i%%4 != 0){
    p2 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
    print(paste(i, "p2"))
    amount = amount + 1}
  
  if((i+1)%%4 == 0){
    p3 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
    print(paste(i, "p3"))
    amount = amount + 1}
  
  if(i%%4 == 0){
    p4 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
    print(paste(i, "p4"))
    amount = 0
    
    png(paste(output,i-3,"-",i ,".png", sep=""))
    grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2)
    dev.off()
    
  }
  if(i == length(IDs) & i%%4 != 0){
    if(amount == 1){
      png(paste(output,i ,".png", sep=""))
      dev.off()
    }
    if(amount == 2){
      png(paste(output,i-1,"-",i ,".png", sep=""))
      grid_arrange_shared_legend(p1, p2, ncol = 2, nrow = 2)
      dev.off()  
    }
    else{
      png(paste(output,i-2,"-",i ,".png", sep=""))
      grid_arrange_shared_legend(p1, p2, p3, ncol = 2, nrow = 2)
      dev.off() 
    }
  }
  print(paste(i, " of ", length(IDs), " done ", sep = ""))
}


#### PCOA of all data ####

phenos <- select(PCOA_df, sample_name, Phenotype, Study_ID)
phenos$Phenotype[phenos$Phenotype=="1"] <- "sink"
phenos$Phenotype[phenos$Phenotype=="0"] <- "non sink"


phenos$Study_ID <- Studies$Author[match(phenos$Study_ID, Studies$ID)] 


phenos$Phenotype <- paste(phenos$Study_ID, phenos$Phenotype)

row.names(PCOA_df) <- PCOA_df$sample_name
numeric_df <- PCOA_df[,6:length(PCOA_df)]
row.names(numeric_df) <- phenos$sample_name



bray <- vegdist(numeric_df, method = "bray")
pcoa_val <- pco(bray, negvals = "zero", dround = 0)

# pco.labels = lapply(row.names(pcoa_val$vectors), function(x) unlist(strsplit(x, split = ".", fixed=TRUE))[1]) # remove .numbers from the end of the row names
# pco.labels = row.names(pcoa_val$vectors) # remove .numbers from the end of the row names

eigenvalues = pcoa_val$values


pcoa_val.df = data.frame(sample_name = row.names(pcoa_val$vectors),
                         PCoA1 = pcoa_val$vectors[,1],
                         PCoA2 = pcoa_val$vectors[,2])

pcoa_val.df = merge(pcoa_val.df, phenos, by.all = "sample_name")

pco.plot = ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(col = as.factor(Phenotype)), alpha = 0.7) +
  scale_color_brewer(name = "phenotype", palette = "Accent") +
  labs(title = "PCoA of Count data", x = paste("PC1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
       y = paste("PC2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""))

png(paste(output,"_lcombine.png",sep=""))
print(pco.plot)
dev.off()

print("script complete")
