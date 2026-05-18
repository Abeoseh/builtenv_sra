.libPaths( c( "~/my_R_libs", .libPaths()) )

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
do_permanovas = args[5]

print(paste("pheno1:", pheno1, "pheno2:", pheno2))
output = strsplit(output, ".png")[[1]]

PCOA_df <- read.csv(input)

Studies = read.csv("./csv_files/phenotypes.csv")


## make an ontology column
meta_files <- list.files(".", "*_SraRunTable.txt", recursive = TRUE)[2:5]

PCOA_df$Ontology <- NA
ontology_column_names <- c("sample_type", "surface_sampledsample", "surfaceaggregated", "surface")
for (i in 1:length(meta_files)){
  ontology_df <- read.csv(meta_files[i])
  ontology_column <- ontology_column_names[i]
  
  matched_ontology_rows <- match(PCOA_df$sample_name, ontology_df$Run)[!is.na(match(PCOA_df$sample_name, ontology_df$Run))]
  matched_pcoa_df_rows <- match(ontology_df$Run, PCOA_df$sample_name)[!is.na(match(ontology_df$Run, PCOA_df$sample_name))]
  
  PCOA_df$Ontology[matched_pcoa_df_rows] <- ontology_df[[ontology_column]][matched_ontology_rows]
  
  
}

# write.csv(PCOA_df, "./PCOA_df.csv", row.names = F)

#### PERMANOVAs ####

#print("Starting PERMANOVAs")



PCOA_df = merge(PCOA_df, Studies, by.x = "Study_ID", by.y = "ID", all.x = T)
PCOA_df = relocate(PCOA_df, Author, Ontology)
colnames(PCOA_df)[1:8]

IDs = unique(PCOA_df$Author)
print("IDs:")
print(IDs)

data_cols <- filter(PCOA_df, Author == IDs[1])
data_cols <- data_cols[,-c(1:6)]

dim(data_cols)
factors <- PCOA_df[,c(1,2,3,4,5,8)] ## Author, Ontology, Study_ID, Phenotype, timepoint



dim(factors)
print(paste("colnames factors:", colnames(factors)))


perm_df = data.frame()

if (!is.na(do_permanovas)){ ## id do_permanovas is not NA then make the permanovas
  print("starting PERMANOVAs")
  
  ## PCOA for each timepoint individually ##
  for(ID in IDs){
  
    # data_cols <- filter(PCOA_df, Author == !!ID)
    # data_cols <- data_cols[,-c(1:8)]
  
    meta = filter(factors, Author == ID)
    timepoints <- unique(meta$timepoint)
  
    for (current_timepoint in timepoints){
      current_timepoint_df_meta = filter(meta, timepoint == current_timepoint)
      current_timepoint_df_data_cols = filter(PCOA_df, Author == !!ID & timepoint == current_timepoint)
      current_timepoint_df_data_cols <- current_timepoint_df_data_cols[,-c(1:8)]
      
      perm_df2 <- adonis2(current_timepoint_df_data_cols ~ Phenotype, current_timepoint_df_meta, permutations = 1000) %>% as.data.frame()
  
      perm_df2$Author = ID
      perm_df2$timepoint = current_timepoint
      perm_df <- rbind(perm_df, perm_df2)
  
    }
    print(paste("done with", ID))
  }
  
  write.csv(perm_df, paste(output,"_individual_timepoint_permanova.csv",sep=""))

  perm_df <- data.frame()

  ## PERMANOVA for phenotype, ontology, and timepoint ##
  for(ID in IDs){
    data_cols <- filter(PCOA_df, Author == !!ID)
    data_cols <- data_cols[,-c(1:8)]
  
    meta = filter(factors, Author == ID)
  
    perm_df2 <- adonis2(data_cols ~ Phenotype, meta, permutations = 1000) %>% as.data.frame()
    perm_df2$factor <- row.names(perm_df2)
    perm_df2$Author = ID
  
    perm_df3 <- adonis2(data_cols ~ Ontology, meta, permutations = 1000) %>% as.data.frame()
    perm_df3$factor <- row.names(perm_df3)
    perm_df3$Author = ID
  
    perm_df4 <- adonis2(data_cols ~ timepoint, meta, permutations = 1000) %>% as.data.frame()
    perm_df4$factor <- row.names(perm_df4)
    perm_df4$Author = ID
  
  
    perm_df <- rbind(perm_df, perm_df2, perm_df3, perm_df4)
  
    print(paste("done with", ID))
  }
  
  write.csv(perm_df, paste(output,"_individual_permanova.csv",sep=""))
  
  ##### PERMANOVA with studies as a random term and phenotype as a fixed term for each timepoint in the hospital study alongside 0013/07
  
  perm_df <- data.frame()
  ## permanova between 0013/07 and the other timepoints
  PRJEB14474_data_cols <- filter(PCOA_df, Study_ID == "PRJEB14474")
  for(timepoint in unique(PRJEB14474_data_cols$timepoint)){
    if (timepoint != "0013/07"){
      data_cols = filter(PRJEB14474_data_cols, timepoint %in% c("0013/07", !!timepoint))
      data_cols <- data_cols[,-c(1:8)]
      
      meta = filter(factors, timepoint %in% c("0013/07", !!timepoint))
      
      perm_df1 <- adonis2(data_cols ~ timepoint, meta, permutations = 1000, strata = meta$Phenotype) %>% as.data.frame()
      perm_df1$factor <- timepoint
      perm_df1$Author = "PRJEB14474"
      perm_df <- rbind(perm_df, perm_df1)
    }
  }
    write.csv(perm_df, paste(output,"_0013_07comparison_permanova.csv",sep=""))
    
  
}

#### make bar chart of R2 ####

perm_df <- read.csv(paste(output,"_individual_permanova.csv",sep=""), check.names = FALSE, row.names = 1)

Authors <- Studies$Author
perm_df <- as.data.frame(perm_df)
for (ID in Authors){

  bar_df <- filter(perm_df, Author == !!ID, !if_any(everything(), is.na))


  bar_df$factor[bar_df$factor == "Phenotype"] <- "Sample Name"
  
  png(paste(output,"_",ID,"bar chart.png"))
  bar <- ggplot(bar_df, mapping = aes(x = R2, y = factor)) + geom_col() + 
    labs(title = paste("R2 chart for", ID), x = "R2", y = "Factor", 
         caption = paste("p-value Sample Name:", signif(as.numeric(bar_df[bar_df$factor == "Sample Name",][["Pr(>F)"]]), 3), 
         "   p-value Ontology:", signif(as.numeric(bar_df[bar_df$factor == "Ontology",][["Pr(>F)"]]), 3),
         "   p-value timepoint:", signif(as.numeric(bar_df[bar_df$factor == "timepoint",][["Pr(>F)"]]), 3), sep="")
         )
  print(bar)
  dev.off()
}




#### POCAs ####
IDs = unique(PCOA_df$Study_ID)

group.colors <- c(setNames("#880808",pheno1), setNames("#333BFF", pheno2))


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
    labs(title = str_wrap(chosen_title,30), 
         x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
         y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), 
         color="Sample Name",
         caption = paste("PERMANOVA p-value =", signif(perm_study[["Pr(>F)"]],2), "R2 =", round(perm_study$R2, 2)) )
  return(pco.plot)
}

#### PCOA for each timepoint colored by phenotype (sample name) ####
print("#### PCOA for each timepoint colored by phenotype (sample name) ####")

perm_df <- read.csv(paste(output,"_individual_timepoint_permanova.csv",sep=""), row.names=1, check.names=F)
num_of_images = 0

for(i in 1:length(IDs)){
  
  post_DEBIAS <- filter(PCOA_df, Study_ID == IDs[i])
  timepoints <- unique(post_DEBIAS$timepoint)
  print(paste("number of timepoints:", length(timepoints)))
  for ( j in 1:length(timepoints) ){

    num_of_images = num_of_images + 1
    post_DEBIAS <- filter(PCOA_df, timepoint == timepoints[j])
    
    Study = Studies$Author[Studies$ID == IDs[i]]
    timepoint = unique(post_DEBIAS$timepoint)

    phenos <- select(post_DEBIAS, sample_name, Phenotype, timepoint)
    phenos$Phenotype[phenos$Phenotype=="1"] <- pheno1
    phenos$Phenotype[phenos$Phenotype=="0"] <- pheno2
    
    # print(colnames(post_DEBIAS)[1:10])
    row.names(post_DEBIAS) <- post_DEBIAS$sample_name
    post_DEBIAS <- post_DEBIAS[,9:length(post_DEBIAS)]
    
    row.names(post_DEBIAS) <- phenos$sample_name

    ## for adding the permanovas to the captions
    perm_study <- filter(perm_df, Author == !!Study & timepoint == !!timepoint & "Phenotype" %in% row.names(perm_df))
    perm_study <- perm_study[1,]
    # print(perm_study)
    # print(Study)
    
    # ex i==1, 1+1%%2 == 0 and 1+1%%4 is not 0
    if((j+1)%%2 == 0 & (j+1)%%4 != 0){
      p1 = pcoa(post_DEBIAS, paste("PCoA on ", Study, " at timepoint: ", timepoint, sep = ""))
      print(paste(j, "p1"))
      amount = 1}
    
    if(j%%2 == 0 & j%%4 != 0){
      p2 = pcoa(post_DEBIAS, paste("PCoA on ", Study, "at timepoint:", timepoint, sep = ""))
      print(paste(j, "p2"))
      amount = amount + 1}
    
    if((j+1)%%4 == 0){
      p3 = pcoa(post_DEBIAS, paste("PCoA on ", Study, "at timepoint:", timepoint, sep = ""))
      print(paste(j, "p3"))
      amount = amount + 1}
    
    if(j%%4 == 0){
      p4 = pcoa(post_DEBIAS, paste("PCoA on ", Study, "at timepoint:", timepoint, sep = ""))
      print(paste(j, "p4"))
      amount = 0
      
      png(paste(output,j-3,"-",j ,"_",num_of_images,"_",Study,"_timepoint_stratified.png", sep=""))
      grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2)
      dev.off()
      
    }
    if(j == length(timepoints) & j%%4 != 0){
      
      if(amount == 1){
        png(paste(output,j ,"_",num_of_images,"_",Study,"_timepoint_stratified.png", sep=""))
        print(p1)
        dev.off()
        
      } else if(amount == 2){
        png(paste(output,j-1,"-",j ,"_",num_of_images,"_",Study,"_timepoint_stratified.png", sep=""))
        grid_arrange_shared_legend(p1, p2, ncol = 2, nrow = 2)
        dev.off() 
        
      } else{
        png(paste(output,j-2,"_",num_of_images,"_",Study,"-",j ,"_timepoint_stratified.png", sep=""))
        grid_arrange_shared_legend(p1, p2, p3, ncol = 2, nrow = 2)
        
        dev.off() 
      }
    }
    
    
  }
  print(paste(i, " of ", length(IDs), " done ", sep = ""))
  
}


#### PCOA for each study colored by all time points ####
print("#### PCOA for each study colored by all time points ####")

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

  ## phenos is made below and will give me the timepoint column
  pcoa_val.df = merge(pcoa_val.df, phenos, by.all = "sample_name")
  pcoa_val.df$timepoint <- as.factor(pcoa_val.df$timepoint)
  
  
  # compute variance explained by each PCoA via eigenvalues 
  eigenvalues = pcoa_val$values
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = timepoint), alpha = 0.7) +
    scale_color_brewer(name = "Timepoint", palette = "Paired") +
    stat_ellipse(level = 0.95, aes(group = timepoint, color = timepoint)) +
    theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
    labs(title = chosen_title, 
         x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
         y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), 
         caption = paste("PERMANOVA p-value =", signif(perm_study[["Pr(>F)"]],2), "R2 =", round(perm_study$R2, 2)) )
  return(pco.plot)
}





perm_df <- read.csv(paste(output,"_individual_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- filter(perm_df, factor == "timepoint")


for(i in 1:length(IDs)){
  post_DEBIAS <- filter(PCOA_df, Study_ID == IDs[i])

  
  phenos <- select(post_DEBIAS, sample_name, timepoint)
  
  row.names(post_DEBIAS) <- post_DEBIAS$sample_name
  post_DEBIAS <- post_DEBIAS[,8:length(post_DEBIAS)]
  row.names(post_DEBIAS) <- phenos$sample_name
  
  
  Study = Studies$Author[Studies$ID == IDs[i]]
  ## for adding the permanovas to the captions
  perm_study <- filter(perm_df, Author == !!Study & "timepoint" %in% row.names(perm_df)) %>% select(-c(factor)) 
  perm_study <- perm_study[1,]
  # print(perm_study)
  # print(Study)
  
  p1 = pcoa(post_DEBIAS[2:length(post_DEBIAS)], paste("PCoA for Timepoint for", Study, sep = ""))
  png(paste(output,"_",IDs[i],"_timepoint_PCOA.png", sep=""))
  print(p1)
  dev.off()

  print(paste(i, " of ", length(IDs), " done ", sep = ""))
}





#### PCOA for all studies together colored by sample name ####
print("#### PCOA for all studies together colored by sample name (Phenotype) ####")

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
  
  group.colors <- c(setNames("#880808",pheno1), setNames("#333BFF", pheno2))

  
  pcoa_val.df = data.frame(sample_name = row.names(pcoa_val$vectors),
                           PCoA1 = pcoa_val$vectors[,1], 
                           PCoA2 = pcoa_val$vectors[,2])
  
  ## phenos is made below and will give me the sample name column (Phenotype)
  pcoa_val.df = merge(pcoa_val.df, phenos, by.all = "sample_name")
  pcoa_val.df$Phenotype <- as.factor(pcoa_val.df$Phenotype)
  
  
  # compute variance explained by each PCoA via eigenvalues 
  eigenvalues = pcoa_val$values
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = Phenotype), alpha = 0.7) +
    
    # scale_color_brewer(name = "Sample Name", palette = "Paired") +
    scale_color_manual(name = "Sample Name", values = group.colors) +
    stat_ellipse(level = 0.95, aes(group = Phenotype, color = Phenotype)) +
    theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
    labs(title = chosen_title, 
         x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
         y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), 
         caption = paste("PERMANOVA p-value =", signif(perm_study[["Pr(>F)"]],2), "R2 =", round(perm_study$R2, 2)) )
  return(pco.plot)
}





perm_df <- read.csv(paste(output,"_individual_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- filter(perm_df, factor == "Phenotype")


for(i in 1:length(IDs)){
  post_DEBIAS <- filter(PCOA_df, Study_ID == IDs[i])
  
  
  phenos <- select(post_DEBIAS, sample_name, Phenotype)
  phenos$Phenotype[phenos$Phenotype=="1"] <- pheno1
  phenos$Phenotype[phenos$Phenotype=="0"] <- pheno2
  
  row.names(post_DEBIAS) <- post_DEBIAS$sample_name
  post_DEBIAS <- post_DEBIAS[,8:length(post_DEBIAS)]
  # row.names(post_DEBIAS) <- phenos$sample_name
  
  
  Study = Studies$Author[Studies$ID == IDs[i]]
  ## for adding the permanovas to the captions
  
  perm_study <- filter(perm_df, Author == !!Study & "Phenotype" %in% row.names(perm_df)) %>% select(-c(factor)) 
  perm_study <- perm_study[1,]
  print(perm_study)
  print(Study)
  
  p1 = pcoa(post_DEBIAS[2:length(post_DEBIAS)], paste("PCoA of Sample Type for ", Study, sep = ""))
  
  png(paste(output,IDs[i],"-","Phenotype_PCOA.png", sep=""))
  print(p1)
  dev.off()
  
  
  print(paste(i, " of ", length(IDs), " done ", sep = ""))
}









#### PCOA between 0013/07 and other timepoints colored by phenotype (sample name) and timepoint ####
print("____________________________________________________________________________________________________________")
print("#### PCOA between 0013/07 and other timepoints colored by phenotype (sample name) and timepoint ####")
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
  pcoa_val.df$Phenotype <- as.factor(pcoa_val.df$Phenotype)
  pcoa_val.df$TimepointPhenotype <- as.factor(pcoa_val.df$TimepointPhenotype)
  print(colnames(pcoa_val.df))
  
  timepoint1 <- unique(pcoa_val.df$timepoint)[unique(pcoa_val.df$timepoint) != "0013/07"]
  # group.colors <- c(setNames("#880808",timepoint1), setNames("#333BFF", "0013/07"))
  group.colors <- c(setNames("#880808",pheno1), setNames("#333BFF", pheno2))
  
  # compute variance explained by each PCoA via eigenvalues 
  eigenvalues = pcoa_val$values
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = Phenotype, shape = timepoint), alpha = 0.7) +
    
    stat_ellipse(level = 0.95, aes(group = TimepointPhenotype, linetype = TimepointPhenotype, )) +
    scale_linetype_manual(name = "Ellipse", values = c(1,2,9,6)) +
    scale_color_manual(values = group.colors) +
    # scale_color_brewer(name = "timepoint - Sample Name", palette = "Paired") +
    
    scale_shape_manual(values=c(3, 16))+
    theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
    labs(title = str_wrap(chosen_title,30), 
         x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
         y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), 
         caption = paste("PERMANOVA p-value =", signif(perm_study[["Pr(>F)"]],2), "R2 =", round(perm_study$R2, 2)) )

  return(pco.plot)
}


perm_df <- read.csv(paste(output,"_0013_07comparison_permanova.csv",sep=""), row.names=1, check.names=F)

PRJEB14474post_DEBIAS <- filter(PCOA_df, Study_ID == "PRJEB14474")
timepoints <- unique(PRJEB14474post_DEBIAS$timepoint)
for(i in 1:length(timepoints)){
  
  if (timepoints[i] != "0013/07"){
    phenos <- select(PRJEB14474post_DEBIAS, sample_name, Phenotype, timepoint) %>% filter(timepoint %in% c("0013/07", timepoints[i]))
    phenos$Phenotype[phenos$Phenotype=="1"] <- pheno1
    phenos$Phenotype[phenos$Phenotype=="0"] <- pheno2
    # phenos$Phenotype <- paste(phenos$timepoint, "-" , phenos$Phenotype)
    phenos$TimepointPhenotype <- paste(phenos$timepoint, "-" , phenos$Phenotype)
    
    ## For naming the graph
    Study = Studies$Author[Studies$ID == "PRJEB14474"]
    
    ## for adding the permanovas to the captions

    perm_study <- filter(perm_df, factor == !!timepoints[i] & "timepoint" %in% row.names(perm_df))
    
    
    perm_study <- perm_study[1,]
    # print(perm_study)
    # print(perm_study[["Pr(>F)"]])
    # print
    # print(Study)
    
    row.names(PRJEB14474post_DEBIAS) <- PRJEB14474post_DEBIAS$sample_name
    post_DEBIAS <- PRJEB14474post_DEBIAS[,9:length(post_DEBIAS)]
    
    p1 = pcoa(post_DEBIAS, paste("PCoA on ", Study, " at timepoint: ", timepoints[i], sep = ""))
    png(paste(output,str_replace(timepoints[i], "/", "-"),"_against_0013-07.png", sep=""))
    print(p1)
    dev.off()
    }
  }




print(warnings())
print("script complete")
