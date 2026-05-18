.libPaths( c( "~/my_R_libs", .libPaths()) )

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




plot_data_columns <- colnames(plot_data)

# print(plot_data_columns)

# print(dim(plot_data))

## sum the bacterial counts
plot_data <- aggregate(plot_data$count, list(plot_data$sample_name, plot_data$Phenotype, plot_data$Study_ID,  
                                             plot_data$timepoint, plot_data$bacteria), FUN=sum)

colnames(plot_data) <- plot_data_columns


## display the top 10 bacteria in each timepoint then leave the rest as other
# ggplot_data <- data.frame()
# for (timepoint in unique(plot_data$timepoint)){
#   
#   top_10 <- plot_data %>% filter(timepoint == !!timepoint) %>% arrange(desc(count))
#   not_top_10_count <- top_10$count[11:length(top_10)] %>% sum()
#   
#   top_10 <- top_10[1:10,] 
#   ggplot_data <- rbind(ggplot_data, top_10, data.frame(timepoint = timepoint, bacteria = "other", count = not_top_10_count))
#   
# }
# 
# ggplot_data$Study_ID <- lognorm_data$Study_ID[match(ggplot_data$timepoint, lognorm_data$timepoint)]
# ggplot_data <- ggplot_data %>% merge(phenotypes, by.x = "Study_ID", by.y = "ID")


### get the top 10 bacteria
top_10 <- aggregate(plot_data$count, list(plot_data$bacteria), FUN=sum) %>% arrange(desc(x)) ## aggregate renames bacteria to Group.1 and count to x

# print(head(top_10))
not_top_10_count <- top_10$x[11:length(top_10)] %>% sum()
top_10 <- top_10$Group.1[1:10]
plot_data$bacteria[!plot_data$bacteria %in% top_10 ] <- "Other"



# print(head(plot_data))


plot_data <- plot_data %>% merge(phenotypes, by.x = "Study_ID", by.y = "ID")
plot_data$Author <- sapply(str_split(plot_data$Author, ":"), `[`, 1) 

# year <- sapply(str_split(plot_data$timepoint, "/"), `[`, 1)
# month <- sapply(str_split(plot_data$timepoint, "/"), `[`, 2)
# plot_data$timepoint <- paste(year, month, sep=" ")
# print(head(plot_data))


## plot the data, generates one graph for each Author

for (Author in unique(plot_data$Author)){
  ggplot_data <- filter(plot_data, Author == !!Author)
  # data$sample_name <- sapply(str_split(data$timepoint, "-"), `[`, 3)
  # data$Phenotype <- sapply(str_split(data$timepoint, "-"), `[`, 2)
  # data$timepoint <- sapply(str_split(data$timepoint, "-"), `[`, 1)
  
  # Define the number of colors you want
  # nb.cols <- dim(data)[1]
  # mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols) ## extend the Dark2 color palette
  
  
  png(paste(output, Author, ".png", sep=""), width = 2000)
  # a <- interaction(ggplot_data$sample_name, ggplot_data$timepoint, ggplot_data$Phenotype, sep = ".")
  # print(a)
  
  plot <- ggplot(ggplot_data) +
  geom_bar(aes(x = interaction(sample_name, timepoint, Phenotype, sep = "."), y = count, fill = bacteria, color = bacteria), stat = "identity", position = "fill") +
  # geom_bar(aes(x = interaction(sample_name, Author), y = count))+
  # scale_x_discrete(NULL, guide = "axis_nested") +
  # guides(x = "axis_nested") +
  guides(x = guide_axis_nested(key_range_auto("\\."))) +
  # scale_fill_manual(values = mycolors) +
  labs(title = paste("Relative Abundance for",Author), x = "Sample Type") +
  theme(legend.text = element_text(size=9), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 10))
  # theme(legend.text = element_text(size=9), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 10),
        # axis.text.x = element_text(angle = 90, hjust = 1)) #+
  # guides(fill = guide_legend(ncol = 1))
  
  print(plot)
  dev.off()
}

print("script complete")
