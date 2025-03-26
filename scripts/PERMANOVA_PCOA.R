.libPaths( c( .libPaths(), "~/my_R_libs") )

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
# PCoA
suppressPackageStartupMessages(library(ecodist))
suppressPackageStartupMessages(library(vegan))
# PCA
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(lemon))
suppressPackageStartupMessages(library(ggplot2))

### IF statement explaination:
# the images are arranged in 2x2 grids

args = commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]
pheno1 =  args[3]
pheno2 = args[4]

output = strsplit(output, ".png")[[1]]

PCOA_df <- read.csv(input)
Studies = read.csv("./csv_files/phenotypes.csv")


#### PERMANOVA ####

print("Starting PERMANOVAs")


PCOA_df = merge(PCOA_df, Studies, by.x = "Study_ID", by.y = "ID", all.x = T)

PCOA_df = relocate(PCOA_df, Author)
colnames(PCOA_df)[1:6]

IDs = unique(PCOA_df$Author)
print("IDs:")
print(IDs)

data_cols <- filter(PCOA_df, Author == IDs[1])
data_cols <- data_cols[,-c(1:6)]


factors <- PCOA_df[,c(1,2,4)]


meta = filter(factors, Author == IDs[1])

perm_df <- adonis2(data_cols ~ Phenotype, meta, permutations = 1000) %>% as.data.frame()
perm_df$Author = IDs[1]
print(paste("done with", IDs[1]))
IDs = IDs[-c(1)]


for(ID in IDs){
  data_cols <- filter(PCOA_df, Author == !!ID)
  data_cols <- data_cols[,-c(1:6)]

  meta = filter(factors, Author == ID)

  perm_df2 <- adonis2(data_cols ~ Phenotype, meta, permutations = 1000) %>% as.data.frame()
  perm_df2$Author = ID
  perm_df <- rbind(perm_df, perm_df2)
  print(paste("done with", ID))
}


write.csv(perm_df, paste(output,"_individual_permanova.csv",sep=""))

#### PERMANOVA with studies as a random term and phenotype as a fixed term


data_cols <- PCOA_df[,-c(1:6)]
print("TRUE %in% is.na(data_cols) ")
print(TRUE %in% is.na(data_cols))
meta <- PCOA_df[,c(1,2,4)]

# perm = adonis2(data_cols ~ Phenotype, method = "bray", meta, permutations = 1000, strata = meta$Study_ID)
# print(perm)
# print(summary(perm))

perm_df <- adonis2(data_cols ~ Phenotype, method = "bray", meta, permutations = 1000, strata = meta$Study_ID) %>% as.data.frame()


write.csv(perm_df, paste(output,"_stratified_permanova.csv",sep=""))


#### PERMANOVA of Study_ID


data_cols <- PCOA_df[,-c(1:7)]
meta <- PCOA_df[,c(1,2,4)]

# perm = adonis2(data_cols ~ Study_ID, method = "bray", meta, permutations = 1000)
# print(perm)
# print(summary(perm))

perm_df <- adonis2(data_cols ~ Study_ID, method = "bray", meta, permutations = 1000) %>% as.data.frame()


write.csv(perm_df, paste(output,"_StudyID_permanova.csv",sep=""))


print("PERMANOVA finished, starting POCAs")

#### POCA ####
IDs = unique(PCOA_df$Study_ID)

group.colors <- c(setNames("#880808",pheno1), setNames("#333BFF", pheno2))

#### PCOA individual ####

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
  pcoa_val.df$Phenotype <- as.factor(pcoa_val.df$Phenotype)

  
  # compute variance explained by each PCoA via eigenvalues 
  eigenvalues = pcoa_val$values
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = Phenotype), alpha = 0.7) +
    scale_color_manual(values = group.colors) +
    stat_ellipse(level = 0.95, aes(group = Phenotype, color = Phenotype)) +
	  theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
    labs(title = str_wrap(chosen_title,30), x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
        y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), color="Phenotype",
        caption = paste("PERMANOVA p-value =", signif(perm_study[["Pr(>F)"]],2), "R2 =", round(perm_study$R2, 2)) )
  return(pco.plot)
}

#### Individual PCOAs ####
perm_df <- read.csv(paste(output,"_individual_permanova.csv",sep=""), row.names=1, check.names=F)
for(i in 1:length(IDs)){
  post_DEBIAS <- filter(PCOA_df, Study_ID == IDs[i])

  
  phenos <- select(post_DEBIAS, sample_name, Phenotype)
  phenos$Phenotype[phenos$Phenotype=="1"] <- pheno1
  phenos$Phenotype[phenos$Phenotype=="0"] <- pheno2
  
  row.names(post_DEBIAS) <- post_DEBIAS$sample_name
  post_DEBIAS <- post_DEBIAS[,7:length(post_DEBIAS)]
  row.names(post_DEBIAS) <- phenos$sample_name
  
  
  Study = Studies$Author[Studies$ID == IDs[i]]
  perm_study <- filter(perm_df, Author == !!Study & "Phenotype" %in% row.names(perm_df))
  perm_study <- perm_study[1,]
  print(perm_study)
  print(Study)

  # ex i==1, 1+1%%2 == 0 and 1+1%%4 is not 0
  if((i+1)%%2 == 0 & (i+1)%%4 != 0){
    p1 = pcoa(post_DEBIAS[2:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
    print(paste(i, "p1"))
    amount = 1}
  
  if(i%%2 == 0 & i%%4 != 0){
    p2 = pcoa(post_DEBIAS[2:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
    print(paste(i, "p2"))
    amount = amount + 1}
  
  if((i+1)%%4 == 0){
    p3 = pcoa(post_DEBIAS[2:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
    print(paste(i, "p3"))
    amount = amount + 1}
  
  if(i%%4 == 0){
    p4 = pcoa(post_DEBIAS[2:length(post_DEBIAS)], paste("PCoA on ", Study, sep = ""))
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


#### PCOA of all data colored by phenotype ####
perm_df <- read.csv(paste(output,"_stratified_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- perm_df[1,]

phenos <- select(PCOA_df, sample_name, Phenotype, Study_ID)
phenos$Phenotype[phenos$Phenotype=="1"] <- pheno1
phenos$Phenotype[phenos$Phenotype=="0"] <- pheno2

## code to merge Study_ID and phenotype column... not needed
# phenos$Study_ID <- Studies$Author[match(phenos$Study_ID, Studies$ID)] 
# phenos$Phenotype <- paste(phenos$Study_ID, phenos$Phenotype)

row.names(PCOA_df) <- PCOA_df$sample_name
numeric_df <- PCOA_df[,7:length(PCOA_df)]
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
  # scale_color_brewer(name = "phenotype", palette = "Accent") +
  stat_ellipse(level = 0.95, aes(group = Phenotype, color = Phenotype)) +
  scale_color_manual(values = group.colors) +
  theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
  labs(title = "PCoA of Count data", x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
       y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""),
       caption = paste("PERMANOVA p-value =", signif(perm_df[["Pr(>F)"]],2), "R2 =", round(perm_df$R2, 2)),
       color="Phenotype" )

png(paste(output,"_combine.png",sep=""))
print(pco.plot)
dev.off()

#### PCOA of all data colored by Study ID ####
perm_df <- read.csv(paste(output,"_StudyID_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- perm_df[1,]

group.colors <- c(`Hospital: Lax et al. 2017` = "#880808", `Air Force: Sharma et al. 2019` = "#333BFF", 
                  `Dorm: Richardson et al. 2019` = "#32a848", `House: Lax et al. 2014` = "#a832a8" ) # #8a8328 burnt yellow color


rownames(PCOA_df) <- make.names(PCOA_df$Study_ID, unique = TRUE)
rownames(PCOA_df) <- gsub("^X", "", rownames(PCOA_df))

numeric_df <- PCOA_df[,7:length(PCOA_df)]
row.names(numeric_df) <- rownames(PCOA_df)

bray <- vegdist(numeric_df, method = "bray")
pcoa_val <- pco(bray, negvals = "zero", dround = 0)

pco.labels = lapply(row.names(pcoa_val$vectors), function(x) unlist(strsplit(x, split = ".", fixed=TRUE))[1]) # remove .numbers from the end of the row names

unique(pco.labels)

eigenvalues = pcoa_val$values


pcoa_val.df = data.frame(Study_ID = unlist(pco.labels),
                         PCoA1 = pcoa_val$vectors[,1],
                         PCoA2 = pcoa_val$vectors[,2])

print(colnames(pcoa_val.df))


pcoa_val.df = merge(pcoa_val.df, Studies, by.x = "Study_ID", by.y = "ID", all.x=TRUE)

print(dim(pcoa_val.df))

pco.plot = ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(col = as.factor(Author)), alpha = 0.7) +
  # scale_color_brewer(name = "phenotype", palette = "Accent") +
  stat_ellipse(level = 0.95, aes(group = Author, color = Author)) +
  scale_color_manual(values = group.colors) +
  theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
  labs(title = "PCoA of Count data", x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
       y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), color="Study",
       caption = paste("PERMANOVA p-value =", signif(perm_df[["Pr(>F)"]],2), "R2 =", round(perm_df$R2, 2)) )

png(paste(output,"_combine_Study_ID.png",sep=""))
print(pco.plot)
dev.off()

print("script complete")
