#!/usr/bin/env/Rscript 

#### Goal: preform a random forest with 17 datasets as testing data and the remaining one as the training data ####
# Do this 101 times and generate 101 ROC curves
# Rscript randomForest.R 3 > run3_out.txt

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
file = args[1]
permutations = as.numeric(args[2])
ID_label = args[3]
folder = args[4]
batch = as.numeric(args[5])

print(ID_label)

DEBIAS_data <- read.csv(file, colClasses = c("Phenotype" = "factor"))

## for naming the variable importance plots
actual_cols <- read.csv(file, colClasses = c("Phenotype" = "factor"), check.names = FALSE)
actual_cols <- names(actual_cols)
colnames(DEBIAS_data)[1:10]

DEBIAS_data %>% distinct(Phenotype) %>% print()

phenos <- read.csv("./csv_files/phenotypes.csv") # for naming the graphs

#### lognorm permutations testing as a loop ####

auc.df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(auc.df) <-   c("Study_ID", "AUC", "Permutation", "DEBIAS")

roc.df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(roc.df) <- c("Study_ID", "DEBIAS", "Permutation", "sensitivities", "specificities")


png(paste("./output/",folder,"/plots/post_DEBIAS-M_RF_lognorm_ROC_", ID_label, ".png", sep=""))#, height = 24, width = 24)


training <- filter(DEBIAS_data, Study_ID != ID_label)
training <- training[, c(2, 4:length(training))]
testing <- filter(DEBIAS_data, Study_ID == ID_label)
testing <- testing[, c(2, 4:length(testing))]

## preform permutations
training_permuted = training
testing_permuted = testing

for(j in 1:permutations){
  ## permutate traing and testing labels

  training_permuted$Phenotype <- sample(training_permuted$Phenotype)

  testing_permuted$Phenotype <- sample(testing_permuted$Phenotype)
  ## do random forest

  RF_fit2 <- randomForest(Phenotype~., data = training_permuted)
  RF_pred <- predict(RF_fit2, testing_permuted, type = "prob")
  rf_roc <- roc(testing_permuted[,1], RF_pred[,1])
  
  ## plot and add data to AUC and ROC dfs
  # if( j == 1){
  #   p <- plot(rf_roc, print.auc=FALSE, add = FALSE)
  # }else{p <- plot(rf_roc, print.auc=FALSE, add = TRUE)}
  
  auc.df[nrow(auc.df) + 1,] = c(ID_label, auc(rf_roc), TRUE, "post")
  
  
  if (batch == 1){df <- data.frame(Study_ID = ID_label, DEBIAS = "post", Permutation = j, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)}else{
    df <- data.frame(Study_ID = ID_label, DEBIAS = "post", Permutation = (j+(permutations*(batch-1))), sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
  }
  
  roc.df <- rbind(roc.df, df)
  
  print(paste(j, " permutations done", sep=""))
}

####
#### predictions of actual data ####
####
## set seed and run rf
if (batch == 1){
  set.seed(100)
  RF_fit <- randomForest(Phenotype~., data = training, importance=TRUE)
  RF_pred <- predict(RF_fit, testing, type = "prob")
  
  ## compute ROC      
  rf_roc <- roc(testing[,1], RF_pred[,1])
  
  ## add AUC from ROC to AUC df
  auc.df[nrow(auc.df) + 1,] = c(ID_label, auc(rf_roc), FALSE, "post")
  auc.df$AUC <- as.numeric(auc.df$AUC) # convert AUC in auc.df to numeric
  
  
  ## add ROC to ROC df
  df <- data.frame(Study_ID = ID_label, DEBIAS = "post", Permutation = 0, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
  roc.df <- rbind(roc.df, df)
  
  
  ## plot ROC of actual data
  # p <- plot(rf_roc, add = TRUE, col = "red", print.auc = TRUE, print.auc.x = 0.17, print.auc.y = 0.09, lwd = 3)
  # phen <- filter(phenos, ID == ID_label)
  # phen <- phen[1,2]
  # title(paste("Training without: ", phen, sep=""), line = + 2.5, cex.main=1.5)
  
  
  
  ## bootstrap p-value 
  # a <- auc.df[auc.df$Permutation == FALSE,]$AUC
  # samp <- auc.df[auc.df$Permutation == TRUE,]$AUC
  # for.pval = length(samp[samp >= a])/1000 
  # 
  # ## add p-value to ROC plot
  # if (for.pval == 0){text(x = 0.1, y = 0.01, label = "p<0.001", col = "red")
  # }else{text(x = 0.1, y = 0.01, label = paste("p=", signif(for.pval, 3),sep=""), col = "red")}
  # 
  # 
  # p
  # 
  # dev.off()

  print("computed AUC of actual data")
}

#### histograms ####

# png(paste("./output/",folder,"/plots/post_DEBIAS-M_RF_lognorm_histogram_", ID_label, ".png", sep=""))
# 
# # a <- auc.df[auc.df$Permutation == FALSE,]$AUC
# # samp <- auc.df[auc.df$Permutation == TRUE,]$AUC
# # z = (a-mean(samp))/sd(samp) # 1 sample z-test
# # for.pval = pnorm(z, lower.tail = FALSE)
# 
# g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
#   geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
#   labs(title = paste("Training without: ", phen , sep=""), y = "count") +
#   annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
#   scale_y_continuous(expand = expansion(mult = c(0, .1)))
# 
# 
# print(g)
# 
# dev.off()

# 
# #### Variable Importance plots ####
# 
# 
# png(paste("./output/",folder,"/plots/post_var_importance_", ID_label, ".png", sep=""), width = 900, height=500)
# par("mar"=c(5,2,5,2))
# randomForest::varImpPlot(RF_fit, 
#                          sort=TRUE, 
#                          main=paste("Variable Importance Plot ", "Training without: ", phen, sep=""))
# 
# dev.off()
# 
# ## as bar chart:
# 
# 
# # make dataframe from importance() output
# actual_cols = actual_cols[-c(1:3)] # Remove Sample_name, Study_ID, and Phenotype
# feat_imp_df <- importance(RF_fit) %>% 
#   data.frame() %>% 
#   mutate(feature = row.names(.))
# 
# feat_imp_df$feature = actual_cols
# 
# 
# feat_imp_df = arrange(feat_imp_df, desc(MeanDecreaseGini))
# feat_imp_df = feat_imp_df[1:50,]
# 
# # plot dataframe
# g <-  ggplot(feat_imp_df, aes(x = reorder(feature, MeanDecreaseGini), 
#                               y = MeanDecreaseGini)) +
#   geom_bar(stat='identity') +
#   coord_flip() +
#   theme_classic() +
#   labs(
#     x     = "Feature",
#     y     = "Importance",
#     title = str_wrap(paste("Variable Importance Plot ", "Training without: ", phen, sep=""),60)
#   )
# 
# png(paste("./output/",folder,"/plots/post_var_importance_bars_", ID_label, ".png", sep=""), width = 680, height=500)
# print(g)
# dev.off()
# 
# ## pval vs inportance plots if pval exists
# if(file.exists( paste("output/", folder, "/pval_v_pval/files/wilcox_pval.csv", sep="") )){
#   suppressPackageStartupMessages(library(ggrepel))
#   
#   wilcox <- read.csv(paste("output/", folder, "/pval_v_pval/files/wilcox_pval.csv", sep=""))
#   wilcox_ID <- paste("pval_", ID_label, sep="")
#   wilcox <- wilcox[, c("bacteria", wilcox_ID)]
#   wilcox <- na.omit(wilcox)
#   
#   wilcox_imp <- merge(wilcox, feat_imp_df, by.x = "bacteria", by.y = "feature", all.x = TRUE)
#   
#   
#   
#   plot <- ggplot(wilcox_imp, aes(x = abs(.data[[wilcox_ID]]), y = .data[["MeanDecreaseGini"]], label = bacteria)) +
#     geom_point() +
#     
#     geom_vline(xintercept=log10(0.05), linetype='dotted', col = 'red') +
# 
#     geom_text_repel(max.overlaps = 10, force_pull = 1, nudge_y = 1,size = 3) +
#     labs(title = paste("log10 p-value vs Feature Importance plot",phen), x = "log10 p-value", y = "Importance") +
#     theme(plot.title = element_text(size=22), axis.text=element_text(size=11),
#           axis.title=element_text(size=15)) 
#   
#   png(paste("./output/",folder,"/plots/pre_pval_imp_", ID_label, ".png", sep=""), width = 1050, height = 480)
#   print(plot)
#   dev.off()
#   
#   
# }

# print(paste(ID_label, " graphs done."))

#### Write to CSV files ####

# CSVs
AUC_filename <- paste("./output/",folder,"/CSVs/builtenv_AUCs.csv",sep="")
# pval_filename <- paste("./output/",folder,"/CSVs/builtenv_AUC_pvals.csv",sep="")
ROC_filename <- paste("./output/",folder,"/CSVs/builtenv_ROCs.csv",sep="")

write.csv(auc.df, paste("./output/",folder,"/CSVs/builtenv_post_DEBIAS_",ID_label,"_batch_",batch,"_AUCs.csv",sep=""), row.names = FALSE)
write.csv(roc.df, paste("./output/",folder,"/CSVs/builtenv_post_DEBIAS_",ID_label,"_batch_",batch,"_ROCs.csv",sep=""), row.names = FALSE)

# ## AUC
# if(file.exists(AUC_filename)){
#   all <- read.csv(AUC_filename)
#   all <- rbind(all, auc.df)
#   write.csv(all, AUC_filename, row.names = FALSE)
# }else{(write.csv(auc.df, AUC_filename, row.names = FALSE))}
# 
# 
# ## ROC
# if(file.exists(ROC_filename)){
#   all <- read.csv(ROC_filename)
#   all <- rbind(all, roc.df)
#   write.csv(all, ROC_filename, row.names = FALSE)
# }else{(write.csv(roc.df, ROC_filename, row.names = FALSE))}



## pvals
# pval.df <- data.frame(Study_ID = ID_label, pval = for.pval, DEBIAS = TRUE)
# 
# if(file.exists(pval_filename)){
#   all <- read.csv(pval_filename)
#   all <- rbind(all, pval.df)
#   write.csv(all, pval_filename, row.names = FALSE)
# }else{(write.csv(pval.df, pval_filename, row.names = FALSE))}


print("done writing to CSVs")
print("done with all post-DEBIAS")


