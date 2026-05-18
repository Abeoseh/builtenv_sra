.libPaths( c( "~/my_R_libs4.5.1", .libPaths()) )

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(legendry))


## Rscript ./scripts/relative_abundance.R "./csv_files/associated_timepoint/lognorm_data.csv" "./output/associated_timepoint/relative_abundance_" "hand associated" floor
## Rscript ./scripts/relative_abundance.R "./csv_files/skinassoc_timepoint/lognorm_data.csv" "./output/skinassoc_timepoint/relative_abundance_" hand "hand associated"
## Rscript ./scripts/relative_abundance.R "./csv_files/skin_floor_timepoint/lognorm_data.csv" "./output/skin_floor_timepoint/relative_abundance_" hand floor




print("started ")
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]
pheno1 <- args[3]
pheno2 <- args[4]
phenotypes <- read.csv("./csv_files/phenotypes.csv")
lognorm_data <- read.csv(input, check.names = F)

# lognorm_data <- read.csv("lognorm_data.csv", check.names = F)
row.names(lognorm_data) <- lognorm_data$sample_name


lognorm_data$Phenotype[lognorm_data$Phenotype=="1"] <- pheno1
lognorm_data$Phenotype[lognorm_data$Phenotype=="0"] <- pheno2



# lognorm_data$timepoint <- paste(lognorm_data$timepoint, "-",lognorm_data$Phenotype, "-" ,lognorm_data$sample_name, sep="")
# 
# plot_data <- lognorm_data %>% select(-c(sample_name, Study_ID, Phenotype, case, ID))  %>% 
#   pivot_longer(!timepoint, names_to = "bacteria", values_to = "count") 


plot_data <- lognorm_data %>% select(-c(case, ID))  %>%
   pivot_longer(-c(sample_name, Phenotype, Study_ID, timepoint), names_to = "bacteria", values_to = "count", )




plot_data <- plot_data %>% merge(phenotypes, by.x = "Study_ID", by.y = "ID")
plot_data$Author <- sapply(str_split(plot_data$Author, ":"), `[`, 1) 
## plot the data, generates one graph for each Author
for (Author in unique(plot_data$Author)){
  
  ggplot_data <- filter(plot_data, Author == !!Author)
  
  ### get the top 20 bacteria
  top_10 <- aggregate(ggplot_data$count, list(ggplot_data$bacteria), FUN=sum) %>% arrange(desc(x)) ## aggregate renames bacteria to Group.1 and count to x
  top_10 <- top_10$Group.1[1:20]
  ggplot_data$bacteria[!ggplot_data$bacteria %in% top_10 ] <- "Other"
  ggplot_data$bacteria <- relevel(as.factor(ggplot_data$bacteria), "Other")
  
  

  
  # ggplot_data$sample_name <- ""  

  default_colors <- scales::hue_pal()(length(unique(ggplot_data$bacteria)))  # ggplot default palette
  names(default_colors) <- unique(ggplot_data$bacteria)
  
  # Override the color for "Other"
  default_colors["Other"] <- "#808080"
  
  png(paste(output, Author, ".png", sep=""), width = 2000)

  plot <- ggplot(ggplot_data) +
    
  geom_bar(aes(x = interaction(sample_name, timepoint, Phenotype, sep = "."), y = count, fill = bacteria, color = bacteria), stat = "identity", position = "fill") +
  ## remove the other samples
  # plot <- ggplot(filter(ggplot_data, bacteria != "Other")) +
  # geom_bar(aes(x = interaction(sample_name, timepoint, Phenotype, sep = "."), y = count, fill = bacteria, color = bacteria), stat = "identity", position = "stack") +
    
  
  labs(title = paste("Relative Abundance for the",Author), x = "Sample Type", y = "") +
    
  theme(legend.text = element_text(size=9), legend.key.size = unit(0.5, 'cm'), 
        legend.title = element_text(size = 10) ) +
          
  guides(x = guide_axis_nested(title = "", type = "box", key = "\\.")) +
    theme(axis.text.x = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 20),
          legend.title=element_text(size=20),
          legend.text = element_text(size=20),
          plot.title = element_text(size=20)) +
  scale_color_manual(values = default_colors) +
  scale_fill_manual(values = default_colors)
  
 
#  guides(x = guide_axis_nested(title = "Sample Type", type = "box", key = "\\.")) 

#  plot <- plot + theme(axis.text.x = element_blank())
  print(plot)
  dev.off()
  
}

print("script complete")
