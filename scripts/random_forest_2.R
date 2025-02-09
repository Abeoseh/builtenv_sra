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


png(paste("./output/",folder,"/ROC_histograms/post_DEBIAS-M_RF_lognorm_ROC_", ID_label, ".png", sep=""))#, height = 24, width = 24)


training <- filter(DEBIAS_data, Study_ID != ID_label)
training <- training[, c(2, 4:length(training))]
testing <- filter(DEBIAS_data, Study_ID == ID_label)
testing <- testing[, c(2, 4:length(testing))]

## set seed and run rf
set.seed(100)
RF_fit <- randomForest(Phenotype~., method = "class", data = training, importance=TRUE)

set.seed(100)
RF_pred <- predict(RF_fit, testing, type = "prob")

## compute ROC
rf_roc <- roc(testing[,1], RF_pred[,1])

## add AUC from ROC to AUC df
auc.df[nrow(auc.df) + 1,] = c(ID_label, auc(rf_roc), FALSE, TRUE)


## add ROC to ROC df
df <- data.frame(Study_ID = ID_label, DEBIAS = TRUE, Permutation = 0, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
roc.df <- rbind(roc.df, df)

## plot ROC of actual data
p <- plot(rf_roc, add = FALSE, col = "red", print.auc = TRUE)
phen <- filter(phenos, ID == ID_label)
phen <- phen[1,2]
title(paste("Training without: ", phen, sep=""), line = + 2.5, cex.main=1.5)

## preform permutations
for(j in 1:permutations){
  ## permutate traing and testing labels
  set.seed(100)
  training$Phenotype <- sample(training$Phenotype)
  set.seed(100)
  testing$Phenotype <- sample(testing$Phenotype)
  set.seed(100)

  ## do random forest
  set.seed(100)
  RF_fit2 <- randomForest(Phenotype~., method = "class", data = training)
  RF_pred <- predict(RF_fit2, testing, type = "prob")
  rf_roc <- roc(testing[,1], RF_pred[,1])

  p <- plot(rf_roc, print.auc=FALSE, add = TRUE)
  auc.df[nrow(auc.df) + 1,] = c(ID_label, auc(rf_roc), TRUE, TRUE)

  df <- data.frame(Study_ID = ID_label, DEBIAS = TRUE, Permutation = j, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
  roc.df <- rbind(roc.df, df)

  print(paste(j, " done post debias", sep=""))
  }
p

dev.off()

# convert AUC in auc.df to numeric
auc.df$AUC <- as.numeric(auc.df$AUC)

png(paste("./output/",folder,"/ROC_histograms/post_DEBIAS-M_RF_lognorm_histogram_", ID_label, ".png", sep=""))

a <- auc.df[auc.df$Permutation == FALSE,]$AUC
samp <- auc.df[auc.df$Permutation == TRUE,]$AUC
z = (a-mean(samp))/sd(samp) # 1 sample z-test
for.pval = pnorm(z, lower.tail = FALSE)

g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
  geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
  labs(title = paste("Training without: ", phen , sep=""), y = "count") +
  annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))


print(g)

dev.off()


#### Variable Importance plots ####


png(paste("./output/",folder,"/ROC_histograms/post_var_importance_", ID_label, ".png", sep=""), width = 900, height=500)
par("mar"=c(5,2,5,2))
randomForest::varImpPlot(RF_fit, 
                         sort=TRUE, 
                         main=paste("Variable Importance Plot ", "Training without: ", phen, sep=""))

dev.off()

## as bar chart:


# make dataframe from importance() output
actual_cols = actual_cols[-c(1:3)] # Remove Sample_name, Study_ID, and Phenotype
feat_imp_df <- importance(RF_fit) %>% 
  data.frame() %>% 
  mutate(feature = row.names(.))

feat_imp_df$feature = actual_cols


feat_imp_df = arrange(feat_imp_df, desc(MeanDecreaseGini))
feat_imp_df = feat_imp_df[1:50,]

# plot dataframe
g <-  ggplot(feat_imp_df, aes(x = reorder(feature, MeanDecreaseGini), 
                              y = MeanDecreaseGini)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_classic() +
  labs(
    x     = "Feature",
    y     = "Importance",
    title = str_wrap(paste("Variable Importance Plot ", "Training without: ", phen, sep=""),60)
  )

png(paste("./output/",folder,"/ROC_histograms/post_var_importance_bars_", ID_label, ".png", sep=""), width = 680, height=500)
print(g)
dev.off()



print(paste(ID_label, " graph done."))

#### Write to CSV files ####

# AUCs
AUC_filename <- paste("./output/",folder,"/AUCs/builtenv_AUCs.csv",sep="")
pval_filename <- paste("./output/",folder,"/AUCs/builtenv_AUC_pvals.csv",sep="")
ROC_filename <- paste("./output/",folder,"/AUCs/builtenv_ROCs.csv",sep="")

write.csv(auc.df, paste("./output/",folder,"/AUCs/builtenv_post_DEBIAS_",ID_label,"_AUCs.csv",sep=""), row.names = FALSE)
write.csv(roc.df, paste("./output/",folder,"/AUCs/builtenv_post_DEBIAS_",ID_label,"_ROCs.csv",sep=""), row.names = FALSE)

## AUC
if(file.exists(AUC_filename)){
  all <- read.csv(AUC_filename)
  all <- rbind(all, auc.df)
  write.csv(all, AUC_filename, row.names = FALSE)
}else{(write.csv(auc.df, AUC_filename, row.names = FALSE))}


## ROC
if(file.exists(ROC_filename)){
  all <- read.csv(ROC_filename)
  all <- rbind(all, roc.df)
  write.csv(all, ROC_filename, row.names = FALSE)
}else{(write.csv(roc.df, ROC_filename, row.names = FALSE))}



## pvals
pval.df <- data.frame(Study_ID = ID_label, pval = for.pval, DEBIAS = TRUE)

if(file.exists(pval_filename)){
  all <- read.csv(pval_filename)
  all <- rbind(all, pval.df)
  write.csv(all, pval_filename, row.names = FALSE)
}else{(write.csv(pval.df, pval_filename, row.names = FALSE))}


print("done writing to CSVs")
print("done with all post-DEBIAS")


