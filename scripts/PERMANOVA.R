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
# pairwise t-test
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# sessionInfo()

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

#### PERMANOVA of Study_ID


data_cols <- PCOA_df[,-c(1:6)]
meta <- PCOA_df[,c(1,2,4)]

perm_df <- adonis2(data_cols ~ Study_ID, method = "bray", meta, permutations = 1000) %>% as.data.frame()
write.csv(perm_df, paste(output,"_StudyID_permanova.csv",sep=""))

## Skin
data <- filter(PCOA_df, Phenotype == 1)
data_cols <- data[,-c(1:6)]
meta <- data[,c(1,2,4)]

perm_df <- adonis2(data_cols ~ Study_ID, method = "bray", meta, permutations = 1000) %>% as.data.frame()
write.csv(perm_df, paste(output,"_StudyID_skin_permanova.csv",sep=""))

## environment
data <- filter(PCOA_df, Phenotype == 0)
data_cols <- data[,-c(1:6)]
meta <- data[,c(1,2,4)]

perm_df <- adonis2(data_cols ~ Study_ID, method = "bray", meta, permutations = 1000) %>% as.data.frame()
write.csv(perm_df, paste(output,"_StudyID_environmental_permanova.csv",sep=""))


print("PERMANOVA finished, starting POCAs")

#### POCA ####
IDs = unique(PCOA_df$Study_ID)

#### PCOA of all data colored by Study ID ####


group.colors <- c(`Hospital: Lax et al. 2017` = "#880808", `Air Force: Sharma et al. 2019` = "#333BFF", 
                  `Dorm: Richardson et al. 2019` = "#32a848", `House: Lax et al. 2014` = "#a832a8" ) # #8a8328 burnt yellow color



## all data 
numeric_df <- PCOA_df[,7:length(PCOA_df)]
row.names(numeric_df) <- PCOA_df$sample_name
# meta = PCOA_df[,c(1:4)]

bray <- vegdist(numeric_df, method = "bray")
pcoa_val <- pco(bray, negvals = "zero", dround = 0)


eigenvalues = pcoa_val$values


pcoa_val.df = data.frame( sample_name = row.names(pcoa_val$vectors),
                         PCoA1 = pcoa_val$vectors[,1],
                         PCoA2 = pcoa_val$vectors[,2] )

print(colnames(pcoa_val.df))


factors = select(PCOA_df, Phenotype, sample_name, Author)
all_pcoa_val.df = merge(pcoa_val.df, factors, by = "sample_name", all.x=TRUE)

## read in the PERMANOVA values
perm_df <- read.csv(paste(output,"_StudyID_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- perm_df[1,]

## plot the PCoA
pco.plot = ggplot(data = all_pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(col = as.factor(Author)), alpha = 0.7) +
  # scale_color_brewer(name = "phenotype", palette = "Accent") +
  stat_ellipse(level = 0.95, aes(group = Author, color = Author)) +
  scale_color_manual(values = group.colors, name = "group") +
  theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
  labs(title = "PCoA of all samples", x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
       y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), color="Study",
       caption = paste("PERMANOVA p-value =", signif(perm_df[["Pr(>F)"]],2), "R2 =", round(perm_df$R2, 2)) )

png(paste(output,"_combine_Study_ID.png",sep=""))
print(pco.plot)
dev.off()


## skin
data <- filter(PCOA_df, Phenotype == 1)
numeric_df <- data[,7:length(data)]
row.names(numeric_df) <- data$sample_name
# meta = PCOA_df[,c(1:4)]

bray <- vegdist(numeric_df, method = "bray")
pcoa_val <- pco(bray, negvals = "zero", dround = 0)


eigenvalues = pcoa_val$values


pcoa_val.df = data.frame( sample_name = row.names(pcoa_val$vectors),
                          PCoA1 = pcoa_val$vectors[,1],
                          PCoA2 = pcoa_val$vectors[,2] )


factors = select(data, Phenotype, sample_name, Author)
pcoa_val_skin.df = merge(pcoa_val.df, factors, by = "sample_name", all.x=TRUE)

perm_df <- read.csv(paste(output,"_StudyID_skin_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- perm_df[1,]

# pcoa_val_skin.df <- filter(pcoa_val.df, Phenotype == 1)
# print("pcoa_val_skin.df Authors")
# distinct(pcoa_val_skin.df, Author) %>% print()
# print(dim(pcoa_val_skin.df))

pco.plot = ggplot(data = pcoa_val_skin.df, mapping = aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(col = as.factor(Author)), alpha = 0.7) +
  stat_ellipse(level = 0.95, aes(group = Author, color = Author)) +
  scale_color_manual(values = group.colors) +
  theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
  labs(title = "PCoA of all hand samples", x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
       y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), color="Study",
       caption = paste("PERMANOVA p-value =", signif(perm_df[["Pr(>F)"]],2), "R2 =", round(perm_df$R2, 2)) )

png(paste(output,"_combine_Study_ID_skin.png",sep=""))
print(pco.plot)
dev.off()


## environmental
data <- filter(PCOA_df, Phenotype == 0)
numeric_df <- data[,7:length(data)]
row.names(numeric_df) <- data$sample_name
# meta = PCOA_df[,c(1:4)]

bray <- vegdist(numeric_df, method = "bray")
pcoa_val <- pco(bray, negvals = "zero", dround = 0)


eigenvalues = pcoa_val$values


pcoa_val.df = data.frame( sample_name = row.names(pcoa_val$vectors),
                          PCoA1 = pcoa_val$vectors[,1],
                          PCoA2 = pcoa_val$vectors[,2] )

factors = select(data, Phenotype, sample_name, Author)
pcoa_val_environmental.df = merge(pcoa_val.df, factors, by = "sample_name", all.x=TRUE)

## read in the PERMANOVA values
perm_df <- read.csv(paste(output,"_StudyID_environmental_permanova.csv",sep=""), row.names=1, check.names=F)
perm_df <- perm_df[1,]

## plot PCoA
pco.plot = ggplot(data = pcoa_val_environmental.df, mapping = aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(col = as.factor(Author)), alpha = 0.7) +
  stat_ellipse(level = 0.95, aes(group = Author, color = Author)) +
  scale_color_manual(values = group.colors) +
  theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(hjust = 0)) +
  labs(title = "PCoA of all environmental samples", x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
       y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""), color="Study",
       caption = paste("PERMANOVA p-value =", signif(perm_df[["Pr(>F)"]],2), "R2 =", round(perm_df$R2, 2)) )

png(paste(output,"_combine_Study_ID_environmental.png",sep=""))
print(pco.plot)
dev.off()

print("PCOAS complete")

#### boxplots ####


print("starting boxplots")

dfs = list(all_pcoa_val.df, pcoa_val_skin.df, pcoa_val_environmental.df)
dfs_names = c("all", "hand", "environmental") 


t_test_plot <- function(pco_df, pco){
  
  # calculate p-values
  pval_df <- stats::pairwise.t.test(pco_df[[pco]], pco_df$Author)$p.value %>% 
    data.frame(check.names = FALSE) %>%
    rownames_to_column("groups") %>%
    pivot_longer(cols = -groups, names_to = "variable", values_to = "p_value") %>% na.omit() %>%
    filter(p_value < 0.05)
  pval_df$p_value <- signif(pval_df$p_value, 3)
  pval_df = mutate(pval_df, groups = purrr::pmap(.l = list(groups, variable), .f = c))
  pval_df$stars <- case_when(pval_df$p_value < 0.05 & pval_df$p_value > 0.01 ~ "*",
                             pval_df$p_value <= 0.01 & pval_df$p_value > 0.001 ~ "**",
                             pval_df$p_value <= 0.001 & pval_df$p_value > 0.0001 ~ "***",
                             pval_df$p_value <= 0.0001 ~ "****")
  
  print(pval_df$stars)
  # print(stats::pairwise.t.test(pco_df[[pco]], pco_df$Author)$p.value )
  
  pco.plot = ggplot(pco_df, aes(Author, .data[[pco]])) +
    labs(title = paste(pco, "for",dfs_names[i], "samples")) +
    geom_boxplot(aes(color = Author)) +
    geom_point(aes(color = Author)) +
    scale_color_manual(values = group.colors, name = "Study") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") +
    if (nrow(pval_df) > 0) {
    geom_signif(comparisons = pval_df$groups, annotations = pval_df$stars, 
                na.rm = TRUE, step_increase = 0.1)}

  
  png(paste(output,"_combine_boxplot_", pco, dfs_names[i], ".png", sep=""))
  print(pco.plot)
  dev.off()
  
  

}

for (i in 1:length(dfs)) {
  
 
  pco1_df <- data.frame(dfs[i]) %>% select(PCoA1, Author) %>% group_by(Author)
  t_test_plot(pco1_df, "PCoA1")

  pco2_df <- data.frame(dfs[i]) %>% select(PCoA2, Author) %>% group_by(Author)
  t_test_plot(pco2_df, "PCoA2")  
  
  # # pairwise t-test
  # # calculate p-values
  # pval_df <- stats::pairwise.t.test(pco1_df$PCoA1, pco1_df$Author)$p.value %>% 
  #   data.frame(check.names = FALSE) %>%
  #   rownames_to_column("groups") %>%
  #   pivot_longer(cols = -groups, names_to = "variable", values_to = "p_value") %>% na.omit() %>%
  #   filter(p_value < 0.05)
  # pval_df$p_value <- signif(pval_df$p_value, 3)
  # pval_df = mutate(pval_df, groups = purrr::pmap(.l = list(groups, variable), .f = c))
  # pval_df$stars <- case_when(pval_df$p_value < 0.05 & pval_df$p_value > 0.01 ~ "*",
  #                            pval_df$p_value <= 0.01 & pval_df$p_value > 0.001 ~ "**",
  #                            pval_df$p_value <= 0.001 & pval_df$p_value > 0.0001 ~ "***",
  #                            pval_df$p_value <= 0.0001 ~ "****")
  # 
  # print(pval_df)
  # print(stats::pairwise.t.test(pco1_df$PCoA1, pco1_df$Author)$p.value )
  # 
  # pco.plot = ggplot(pco1_df, aes(Author, PCoA1, color = Author)) +
  #   labs(title = paste("PCo1 for",dfs_names[i])) +
  #   geom_boxplot() +
  #   geom_point() +
  #   scale_color_manual(values = group.colors) +
  #   geom_signif(comparisons = pval_df$groups, annotations = pval_df$stars, 
  #               na.rm = TRUE, step_increase = 0.4) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # 
  # png(paste(output,"_combine_boxplot_pco1", dfs_names[i], ".png", sep=""))
  # print(pco.plot)
  # dev.off()
  
  
  
  
  
  # print("PCO2")
  # pco2_df <- data.frame(dfs[i]) %>% select(PCoA2, Author) %>% group_by(Author)
  # 
  # # calculate p-values
  # pval_df <- stats::pairwise.t.test(pco1_df$PCoA1, pco1_df$Author)$p.value %>% 
  #   data.frame() %>%
  #   rownames_to_column("groups") %>%
  #   pivot_longer(cols = -groups, names_to = "variable", values_to = "p_value") %>% na.omit() %>%
  #   filter(p_value < 0.05)
  # pval_df$p_value <- signif(pval_df$p_value, 3)
  # pval_df = mutate(pval_df, values = purrr::pmap(.l = list(groups, variable), .f = c))
  # 
  # # plot PCoA
  # pco.plot = ggplot(pco2_df, aes(Author, PCoA2)) +
  #   labs(title = paste("PCo2 for",dfs_names[i])) +
  #   scale_color_manual(values = group.colors) +
  #   geom_boxplot() +
  #   geom_point() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # 
  # png(paste(output,"_combine_boxplot_pco2", dfs_names[i], ".png", sep=""))
  # print(pco.plot)
  # dev.off()
  
  print(paste("done with", dfs_names[i]))
}


print("script complete")


