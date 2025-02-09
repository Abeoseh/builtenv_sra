library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(gridExtra)
library(ecodist)
library(vegan)
library(factoextra)
library(lemon)

#### ensure working directory is correct ####

getwd()
setwd("C:/Users/brean/Downloads/masters/Fodor/builtenv_sra")
getwd()

# input folder, located within ./output
input = "associated"
input = "skin_skinassoc"



#### read in files ####


AUCs <- read.csv(paste("./output/",input,"/AUCs/builtenv_AUCs.csv",sep=""), colClasses = c("DEBIAS" = "factor"))

# post_DEBIAS_files <- list.files("./cluster_runs/folder/DEBIAS-M_runs/",
                          # pattern = "debiased_lognorm", full.names = TRUE) # debiased lognorm for pocas

# debias_weights_files <- list.files(paste("./output/",input,"/DEBIAS-M_runs/",sep=""),
                                   # pattern = "debias_weights", full.names = TRUE)
phenos <- read.csv("./csv_files/phenotypes.csv")


pval.df <- read.csv(paste("./output/",input,"/AUCs/builtenv_AUC_pvals.csv",sep=""), colClasses = c("DEBIAS" = "factor"))

IDs <- distinct(AUCs, Study_ID)$Study_ID


#### Optional: Fix AUC df ####
# if there is a problems with the combined AUCs use this to combine the individual AUCs
a = list.files("./cluster_runs/folder/csv_files/AUCs/semi_leaky_cancer/",
               pattern = "cancer_post_", full.names = TRUE)

post_AUCs = read.csv(a[1])
for(i in 2:length(a)){
  d = read.csv(a[i])
  post_AUCs = rbind(post_AUCs, d)
}

AUCs = rbind(AUCs, post_AUCs)
write.csv(AUCs, "./cluster_runs/folder/csv_files/AUCs/semi_leaky_cancer/cancer_semi_leaky_AUCs.csv", row.names = FALSE)


#### Optional: histograms and p-val df ####
# can make pval.df or histograms

pval.df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(pval.df) = c("Study_ID", "pval", "DEBIAS")

selected_AUCs = filter(AUCs, DEBIAS == 0)

for(i in 1:length(IDs)){
  
  # phen <- phenos[phenos$ID == IDs[i],]$phenotype
  
  auc.df <- filter(selected_AUCs, Study_ID == IDs[i])

  a <- auc.df[auc.df$Permutation == 0,]$AUC
  samp <- auc.df[auc.df$Permutation == 1,]$AUC
  z = (a-mean(samp))/(sd(samp)/sqrt(1))
  for.pval = pnorm(z, lower.tail = FALSE)
  
  pval.df[nrow(pval.df) + 1,] = c(IDs[i], for.pval, 3)
  
  # make histograms
  # png(paste("./cluster_runs/DEBIAS-M_100_perm_fixed/output/pre_DEBIAS-M_RF_lognorm_histogram_fixed_", IDs[i], ".png", sep=""), width = 480, height = 480)
  # g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
  #   geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
  #   labs(title = paste("Training without: ", IDs[i], " (", phen, ")", sep=""), y = "count") +
  #   annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
  #   scale_y_continuous(expand = expansion(mult = c(0, .1)))
  # 
  # 
  # print(g)
  # dev.off()
}

# pval.df <- merge(pval.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
# # pval.df$phenotype[grep("cancer", pval.df$phenotype)] <- "cancer"
# write.csv(pval.df, "./csv_files/folder/AUC_pvals_autism.csv", row.names = FALSE)



#### save histograms and ROCS as a png ####

hist_pdf <- function(input_folder, pattern, output_folder){ 
  # list all the files in the input folder while displaying the full name
  all_images <- list.files(input_folder, full.names = TRUE)
  # get the indices of all file names that match the supplied pattern
  ind <- grep(pattern, all_images)
  # subscript all_images based on that
  all_images <- all_images[ind]

  plots <- lapply(ll <- all_images, function(x){
    img <- as.raster(readPNG(x))
    rasterGrob(img, interpolate = FALSE)
  })
  if (length(all_images) == 2){ggsave(output_folder, marrangeGrob(grobs=plots, nrow=1, ncol=2), width = 4, height = 4, dpi = 300)}
  else {ggsave(output_folder, marrangeGrob(grobs=plots, nrow=2, ncol=2), width = 4, height = 4, dpi = 300)}
  
}

# pre, ROC
hist_pdf(paste("./output/",input,"/ROC_histograms",sep="") ,"pre_DEBIAS-M_RF_lognorm_ROC", paste("./output/",input,"/post_cluster_pngs/ROC_pre_100_perm.png",sep=""))

# pre hist
hist_pdf(paste("./output/",input,"/ROC_histograms",sep=""), "pre_DEBIAS-M_RF_lognorm_hist", paste("./output/",input,"/post_cluster_pngs/hist_pre_100_perm.png",sep=""))

# pre var importance BARS
hist_pdf(paste("./output/",input,"/ROC_histograms",sep=""), "pre_var_importance_bars", paste("./output/",input,"/post_cluster_pngs/pre_var_impor_bars.png",sep=""))

# pre var importance 
hist_pdf(paste("./output/",input,"/ROC_histograms",sep=""), "pre_var_importance_P", paste("./output/",input,"/post_cluster_pngs/pre_var_impor.png",sep=""))

# post, ROC
hist_pdf(paste("./output/",input,"/ROC_histograms",sep=""), "post_DEBIAS-M_RF_lognorm_ROC", paste("./output/",input,"/post_cluster_pngs/ROC_post_100_perm.png",sep=""))

# post hist
hist_pdf(paste("./output/",input,"/ROC_histograms",sep=""), "post_DEBIAS-M_RF_lognorm_hist", paste("./output/",input,"/post_cluster_pngs/hist_post_100_perm.png",sep=""))

# pre var importance BARS
hist_pdf(paste("./output/",input,"/ROC_histograms",sep=""), "post_var_importance_bars", paste("./output/",input,"/post_cluster_pngs/post_var_impor_bars.png",sep=""))

# pre var importance 
hist_pdf(paste("./output/",input,"/ROC_histograms",sep=""), "post_var_importance_P", paste("./output/",input,"/post_cluster_pngs/post_var_impor.png",sep=""))






#### box plots and t-test ####
## Don't modify if you only have 2 samples!! ##

# auc box plot

# t-test between baseline and semi-leaky: t = 0.03496, df = 3, p-value = 0.9743

t_auc = AUCs %>% filter(Permutation == FALSE) %>% select(DEBIAS, AUC) 
with(t_auc, t.test(AUC[DEBIAS == FALSE], AUC[DEBIAS == TRUE], alternative = "two.sided", paired = TRUE))

# next three lines are only needed if i have phenotypes to add...I do
# AUCs$DEBIAS <- as.numeric(levels(AUCs$DEBIAS))[(AUCs$DEBIAS)]
selected_AUCs <- AUCs %>% filter(Permutation == FALSE) %>% merge(y = phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
# selected_AUCs$phenotype[grep("cancer", selected_AUCs$phenotype)] <- "cancer"


auc_b <- ggplot(selected_AUCs, mapping = aes(x=as.factor(DEBIAS), y=AUC)) + 
  geom_boxplot(fill = "darkgrey") + 
  geom_point(aes(color = as.factor(Author)), alpha = 1) +
 scale_color_brewer(name = "Author", palette = "Paired", guide = "none") +
  scale_x_discrete(labels = c("Before DEBIAS-M", "After DEBIAS-M")) +
  labs(x = "", title = "AUCs before and after DEBIAS-M") +
  ylim(0.4, 1.0) +
#  geom_signif(annotation = "*",
#              y_position = 0.95, xmin = 2, xmax = 3) +
  theme_minimal()

auc_b

# p-val box plot


# t-test between baseline and semi-leay: t = 1, df = 3, p-value = 0.391

t_pval = pval.df %>% select(DEBIAS, pval) 
with(t_pval, t.test(pval[DEBIAS == FALSE], pval[DEBIAS == TRUE], alternative = "two.sided", paired = TRUE))

# pval.df$DEBIAS <- as.numeric(levels(pval.df$DEBIAS))[(pval.df$DEBIAS)]
pval.df <- pval.df %>% merge(y = phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE) %>% arrange(desc(DEBIAS))
#pval.df$phenotype[grep("cancer", selected_AUCs$phenotype)] <- "cancer"

pval_b <- ggplot(pval.df, mapping = aes(x = as.factor(DEBIAS), y = log(pval, base = 10))) + 
  geom_boxplot(fill = "darkgrey") + 
  geom_point(aes(color = as.factor(Author)), alpha = 1) +
  scale_color_brewer(name = "Author", palette = "Paired") +
  scale_x_discrete(labels = c("Before DEBIAS-M", "After DEBIAS-M")) +
  labs(x = "", y = "log10 p-value", title = "p-values before and after DEBIAS-M") +
#  geom_signif(annotation = "**",
#              y_position = 1, xmin = 2, xmax = 3) +
  theme_minimal() 

pval_b

png(paste("./output/",input,"/post_cluster_pngs/box_AUC_pval.png",sep=""), width = 1050, height = 480)
grid.arrange(auc_b, pval_b, nrow = 1)
dev.off()


#### pcoa ####

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
  
  pco.labels = lapply(row.names(pcoa_val$vectors), function(x) unlist(strsplit(x, split = ".", fixed=TRUE))[1]) # remove .numbers from the end of the row names
  pcoa_val.df = data.frame(Study_ID = unlist(pco.labels),
                           PCoA1 = pcoa_val$vectors[,1], 
                           PCoA2 = pcoa_val$vectors[,2])
  
  pcoa_val.df = merge(pcoa_val.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
  pcoa_val.df$phenotype[grep("cancer", pcoa_val.df$phenotype)] = "cancer"
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = as.factor(phenotype)), alpha = 0.7) +
    scale_color_brewer(name = "phenotype", palette = "Paired") +
    labs(title = chosen_title, x = "PC1", y = "PC2")
  return(pco.plot)
}




# after debias-m count pcoa

# pd = post_DEBIAS_files

for(i in 1:length(post_DEBIAS_files)){
  post_DEBIAS <- read.csv(post_DEBIAS_files[i])
  rownames(post_DEBIAS) <- make.names(post_DEBIAS$Study_ID, unique = TRUE)
  rownames(post_DEBIAS) <- gsub("^X", "", rownames(post_DEBIAS))
  
  if((i+1)%%2 == 0 & (i+1)%%4 != 0){
    p1 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " post count data.", sep = ""))
    print(paste(i, "p1"))
    amount = 1}

  if(i%%2 == 0 & i%%4 != 0){
    p2 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " post count data", sep = ""))
    print(paste(i, "p2"))
    amount = amount + 1}

  if((i+1)%%4 == 0){
    p3 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " post count data", sep = ""))
    print(paste(i, "p3"))
    amount = amount + 1}
  
  if(i%%4 == 0){
    p4 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " post count data", sep = ""))
    print(paste(i, "p4"))
    amount = 0
    
    png(paste("./output/subset/pcoa_post_DEBIAS",i-3,"-",i ,".png", sep=""))
    grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2)
    dev.off()
    
   }
  if(i == length(post_DEBIAS_files) & i%%4 != 0){
    if(amount == 1){
    png(paste("./output/autism/pcoa_post_DEBIAS",i ,".png", sep=""))
    dev.off()
    }
    if(amount == 2){
      png(paste("./output/autism/pcoa_post_DEBIAS",i-1,"-",i ,".png", sep=""))
      grid_arrange_shared_legend(p1, p2, ncol = 2, nrow = 2)
      dev.off()  
    }
    else{
    png(paste("./output/autism/pcoa_post_DEBIAS",i-2,"-",i ,".png", sep=""))
    grid_arrange_shared_legend(p1, p2, p3, ncol = 2, nrow = 2)
    dev.off() 
    }
  }
  print(paste(i, " of ", length(pd), " done ", sep = ""))
}

#### PCA weights ####

pca_plot <- function(df, chosen_title, first_PC, second_PC){
  # pca
  dw_pca <- prcomp(scaled)
  var = get_pca_var(dw_pca)
  # get the contributions as a df
  dwpcacontr.df = as.data.frame(var$contrib)
  # get percent variation
  pca_contr = summary(dw_pca)$importance[2,]
  
  dwpca.df = as.data.frame(dw_pca$x)
  dwpca.df$Study_ID = row.names(dwpca.df)
  dwpca.df = merge(x = dwpca.df, y = phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
  
  dwpca.plot = ggplot(dwpca.df, aes(x = .data[[first_PC]], y = .data[[second_PC]])) +
    geom_point(aes(col = as.factor(Author))) +
    scale_color_brewer(name = "Author", palette = "Paired") +
    labs(title = chosen_title, x = paste(first_PC, " (", pca_contr[[first_PC]]*100, ") ", sep = ""), y = paste(second_PC, " (", pca_contr[[second_PC]]*100,") ", sep=""))
  
}




for(i in 1:length(debias_weights_files)){
  
  # get the ID
  name = strsplit(debias_weights_files[i], "/")[[1]][5] 
  name = strsplit(name, split = "_")[[1]][4]
  name = strsplit(name, split = "[.]")[[1]][1]  
  name = phenos$Author[phenos$ID == name]
    
  debias_weights <- read.csv(debias_weights_files[i], check.names = FALSE, row.names = 1)
  debias_weights = t(debias_weights)
  debias_weights = data.frame(debias_weights)
  
  scaled <- as.data.frame(scale(debias_weights))
  
  if((i+1)%%2 == 0 & (i+1)%%4 != 0){
    p1 = pca_plot(scaled, paste(name,": PCA on weights", sep = ""), "PC1", "PC2")
    print(paste(i, "p1"))
    amount = 1}
  
  if(i%%2 == 0 & i%%4 != 0){
    p2 = pca_plot(scaled, paste(name,": PCA on weights", sep = ""), "PC1", "PC2")
    print(paste(i, "p2"))
    amount = amount + 1}
  
  if((i+1)%%4 == 0){
    p3 = pca_plot(scaled, paste(name,": PCA on weights", sep = ""), "PC1", "PC2")
    print(paste(i, "p3"))
    amount = amount + 1}
  if(i%%4 == 0){
    p4 = pca_plot(scaled, paste(name,": PCA on weights", sep = ""), "PC1", "PC2")
    print(paste(i, "p4"))
    amount = 0
    
    png(paste("./output/associated_with_2192_renamed/pca_post_DEBIAS_",i-3,"-",i ,".png", sep=""))
    grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2)
    dev.off()
  }  
# writing to files
  if(i == length(debias_weights_files) & i%%4 != 0){
    if(amount == 1){
      png(paste("./output/associated_with_2192_renamed/pca_post_DEBIAS_",i ,".png", sep=""))
      print(p1)
      dev.off()
    }
    if(amount == 2){
      png(paste("./output/associated_with_2192_renamed/pca_post_DEBIAS_",i-1,"-",i ,".png", sep=""))
      grid_arrange_shared_legend(p1, p2, ncol = 2, nrow = 2)
      dev.off()  
    }
    if(amount == 3){
      png(paste("./output/associated_with_2192_renamed/pca_post_DEBIAS_",i-2,"-",i ,".png", sep=""))
      grid_arrange_shared_legend(p1, p2, p3, ncol = 2, nrow = 2)
      dev.off() 
    }}

  print(paste(i, " of ", length(debias_weights_files), " done ", sep = ""))
}


for(i in 1:length(dw)){
  debias_weights = t(debias_weights)
  debias_weights = data.frame(debias_weights)
  # rownames(debias_weights) <- debias_weights$Study_ID
  # debias_weights = debias_weights[,2:length(debias_weights)]
  scaled <- as.data.frame(scale(debias_weights))
  # verify mean == 0 and sd == 1
  # sapply(scaled, mean)
  # sapply(scaled, sd)
  
  
  dw_pca <- prcomp(scaled)
  var = get_pca_var(dw_pca)
  dwpcacontr.df = as.data.frame(var$contrib)
  pca_contr = summary(dw_pca)$importance[2,]
  
  dwpca.df = as.data.frame(dw_pca$x)
  dwpca.df$Study_ID = row.names(dwpca.df)
  dwpca.df = merge(dwpca.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
  dwpca.df$phenotype[grep("cancer", dwpca.df$phenotype)] = "cancer"
  
  dwpca.plot = ggplot(dwpca.df, mapping = aes(x = PC1, y = PC2)) +
    geom_point(aes(col = as.factor(phenotype))) +
    scale_color_brewer(name = "phenotype", palette = "Paired") +
    labs(title = "PCA of weights", x = paste("PC1 (", pca_contr[1], ")", sep = ""), y = paste("PC2 (", pca_contr[2],")", sep=""))
  
  png(paste("./output/semi_leaky/PCA_weights_PC1-2_leaky",IDs[i],".png"))
  print(dwpca.plot)
  dev.off()
  print(paste(i, " of ", length(dw), " done.", sep = ""))
  
}

# same plot but with study labels instead
dwpca.df <- read.csv(debias_weights_files[1], check.names = FALSE, row.names = 1)
debias_weights = t(dwpca.df)
debias_weights = data.frame(debias_weights)
scaled <- as.data.frame(scale(debias_weights))

dw_pca <- prcomp(scaled)
var = get_pca_var(dw_pca)
dwpcacontr.df = as.data.frame(var$contrib)
pca_contr = summary(dw_pca)$importance[2,]

dwpca.df = as.data.frame(dw_pca$x)
dwpca.df$Study_ID = row.names(dwpca.df)
dwpca.df = merge(dwpca.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
# dwpca.df$phenotype[grep("cancer", dwpca.df$phenotype)] = "cancer"

dwpca_labels.plot = ggplot(dwpca.df, mapping = aes(x = PC1, y = PC2)) +
  geom_point(color = "white") +
  geom_text(mapping = aes(label = Study_ID)) +
  labs(title = "PCA of weights")
dwpca_labels.plot

# png("./output/test_permutated/PCA_weights_PC1-2_leaky.png")
# dwpca.plot
# dev.off()



#### make plots of percentage of data missing ####
ggplot(data = missing_no_na) + 
  geom_histogram(mapping = aes(x = percent_missing), bins = 40) +
  labs(title = "percent missing when NAs are coded as 0s")

bacteria = filter(missing_no_na, percent_missing < 0.25)


png("./output/percent_missing.png")
ggplot(data = missing_with_na) + 
  geom_histogram(mapping = aes(x = percent_missing), bins = 40) +
  labs(title = "percent missing")
dev.off()

bacteria1 = filter(missing_with_na, percent_missing < 0.25)

print(FALSE %in% (bacteria1 == bacteria))
