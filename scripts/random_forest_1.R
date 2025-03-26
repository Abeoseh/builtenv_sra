#!/usr/bin/env/Rscript 

#### Goal: preform a random forest with 17 datasets as testing data and the remaining one as the training data ####
# Do this 18 times and generate 18 ROC curves
# Rscript randomForest.R 3 > run3_out.txt

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])
permutations = as.numeric(args[2])
input = args[3]
output = args[4]

lognorm_out <- read.csv(paste("./csv_files/",input,"/lognorm_data.csv",sep=""), colClasses = c("Phenotype" = "factor"))
lognorm_out <- lognorm_out[,c(1:3,6:length(lognorm_out))] # remove case and ID columns

## for naming the variable importance plots
actual_cols <- read.csv(paste("./csv_files/",input,"/lognorm_data.csv",sep=""), colClasses = c("Phenotype" = "factor"), check.names = FALSE)
actual_cols <- actual_cols[,c(1:3,6:length(actual_cols))] # remove case and ID columns
actual_cols <- names(actual_cols)

phenos <- read.csv("./csv_files/phenotypes.csv") # for naming the graphs


#### lognorm permutations testing as a loop ####

IDs <- distinct(lognorm_out, Study_ID)$Study_ID

## create AUC and ROC df
auc.df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(auc.df) <-   c("Study_ID", "AUC", "Permutation", "DEBIAS")


roc.df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(roc.df) <- c("Study_ID", "DEBIAS", "Permutation", "sensitivities", "specificities")



png(paste("./output/",output,"/plots/pre_DEBIAS-M_RF_lognorm_ROC_", IDs[i], ".png", sep=""))#, height = 24, width = 24)
# par(mar=c(3,3,1,0))#, mfrow=c(2,2))


training <- filter(lognorm_out, Study_ID != IDs[i])
training <- training[, c(2, 4:length(training))] # this is okay because case and ID columns are removed.
testing <- filter(lognorm_out, Study_ID == IDs[i])
testing <- testing[, c(2, 4:length(testing))]

## set seed and run rf
set.seed(100)
RF_fit <- randomForest(Phenotype~., method = "class", data = training, importance=TRUE)

## predictions
set.seed(100)    
RF_pred <- predict(RF_fit, testing, type = "prob")

## compute ROC      
rf_roc <- roc(testing[,1], RF_pred[,1])

## add AUC from ROC to AUC df
auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), FALSE, FALSE)

## add ROC to ROC df
df <- data.frame(Study_ID = IDs[i], DEBIAS = FALSE, Permutation = 0, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
roc.df <- rbind(roc.df, df)


## plot ROC of actual data
p <- plot(rf_roc, add = FALSE, col = "red", print.auc = TRUE)
phen <- filter(phenos, ID == IDs[i])
phen <- phen[1,2]
title(paste("Training without: ", phen, sep=""), line = + 2.5, cex.main=1.5)
 
 ## preform permutations     
 for(j in 1:permutations){
     ## permutate traing and testing labels
     set.seed(100)
     training$Phenotype <- sample(training$Phenotype)
     set.seed(100)
     testing$Phenotype <- sample(testing$Phenotype)
     ## do random forest
     set.seed(100)
     RF_fit2 <- randomForest(Phenotype~., method = "class", data = training)
     RF_pred <- predict(RF_fit2, testing, type = "prob")
     rf_roc <- roc(testing[,1], RF_pred[,1])

     ## plot and add data to AUC and ROC dfs
     p <- plot(rf_roc, print.auc=FALSE, add = TRUE)

     auc.df[nrow(auc.df) + 1,] = c(IDs[i], auc(rf_roc), TRUE, FALSE)

     
     df <- data.frame(Study_ID = IDs[i], DEBIAS = FALSE, Permutation = j, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
     roc.df <- rbind(roc.df, df)

     print(paste(j, " done", sep=""))
     }
p

dev.off()  

auc.df$AUC <- as.numeric(auc.df$AUC)

#### Histograms ####
png(paste("./output/",output,"/plots/pre_DEBIAS-M_RF_lognorm_histogram_", IDs[i], ".png", sep=""))

a <- auc.df[auc.df$Permutation == FALSE,]$AUC
samp <- auc.df[auc.df$Permutation == TRUE,]$AUC
z = (a-mean(samp))/(sd(samp)/sqrt(1))
for.pval = pnorm(z, lower.tail = FALSE)

g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
  geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
  labs(title = paste("Training without: ", phen, sep=""), y = "count") +
  annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) 


print("done")  

print(g)

dev.off()  


#### Variable Importance plots ####

# png(paste("./output/",output,"/plots/pre_var_importance_", IDs[i], ".png", sep=""), width = 900, height=500)
# randomForest::varImpPlot(RF_fit, 
#                          sort=TRUE, 
#                          main=paste("Variable Importance Plot ", "Training without: ", phen, sep=""))
# dev.off()

## as bar chart:

# make dataframe from importance() output
actual_cols = actual_cols[-c(1:3)] # Remove Sample, Study_ID, Phenotype

feat_imp_df <- importance(RF_fit) %>% 
  data.frame() %>% 
  mutate(feature = row.names(.)) 


feat_imp_df$feature = actual_cols
feat_imp_df$abs_MeanDecreaseGini <- abs(feat_imp_df$MeanDecreaseGini) 

feat_imp_df = arrange(feat_imp_df, desc(abs_MeanDecreaseGini))


feat_imp_top_50 = feat_imp_df[1:50,]


# plot dataframe
g <-  ggplot(feat_imp_top_50, aes(x = reorder(feature, abs_MeanDecreaseGini), 
                              y = abs_MeanDecreaseGini)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Feature",
    y = "Importance",
    title = str_wrap(paste("Variable Importance Plot ", "Training without: ", phen, sep=""),40)) +
  theme(plot.title = element_text(size=15), axis.text=element_text(size=11),
        axis.title=element_text(size=15)) 

png(paste("./output/",output,"/plots/pre_var_importance_bars_", IDs[i], ".png", sep=""), width = 680, height=600)
print(g)
dev.off()

## pval vs inportance plots if pval exists
if(file.exists( paste("output/", output, "/pval_v_pval/files/wilcox_pval.csv", sep="") )){
  # suppressPackageStartupMessages(library(ggrepel))
  
  wilcox <- read.csv(paste("output/", output, "/pval_v_pval/files/wilcox_pval.csv", sep=""))
  wilcox_ID <- paste("pval_", IDs[i], sep="")
  wilcox <- wilcox[, c("bacteria", wilcox_ID)]
  wilcox <- na.omit(wilcox)
  
  wilcox_imp <- merge(wilcox, feat_imp_df, by.x = "bacteria", by.y = "feature", all.x = TRUE)
  wilcox_imp$MeanDecreaseGini <- log10(wilcox_imp$MeanDecreaseGini + 1)
  
  cor <- cor.test(wilcox_imp$MeanDecreaseGini, wilcox_imp[[wilcox_ID]])
  print(cor)

  
  plot <- ggplot(wilcox_imp, aes(x = abs(.data[[wilcox_ID]]), y = .data[["MeanDecreaseGini"]], label = bacteria)) +
    geom_point() +
    # geom_vline(xintercept=log10(0.05), linetype='dotted', col = 'red') +
    geom_vline(xintercept=-log10(0.05), linetype='dotted', col = 'red') +
    # geom_text_repel(force_pull = 2) +
    annotate("label", x=max(wilcox_imp[[wilcox_ID]]), y=max(wilcox_imp[["MeanDecreaseGini"]])-.01, size = 3, label = paste("p= ", signif(cor$p.value, digits=3), sep="")) +
    annotate("label", x=max(wilcox_imp[[wilcox_ID]]), y=max(wilcox_imp[["MeanDecreaseGini"]])-.05, size = 3, label = paste("Tau= ", signif(cor$estimate,digits=3), sep="")) +
    
    labs(title = paste("log10 p-value vs Feature Importance plot",phen), x = "log10 p-value", y = "log10 Importance") +
    theme(plot.title = element_text(size=22), axis.text=element_text(size=11),
          axis.title=element_text(size=15)) 
  
  png(paste("./output/",output,"/plots/pre_pval_imp_", IDs[i], ".png", sep=""), width = 1050, height = 480)
  print(plot)
  dev.off()
  
  
}



print(paste(i, "of", length(IDs), " graphs done.", IDs[i]))


#### Write to CSV files ####
AUC_filename <- paste("./output/",output,"/CSVs/builtenv_AUCs.csv",sep="")
pval_filename <- paste("./output/",output,"/CSVs/builtenv_AUC_pvals.csv",sep="")
ROC_filename <- paste("./output/",output,"/CSVs/builtenv_ROCs.csv",sep="")


write.csv(auc.df, paste("./output/",output,"/CSVs/builtenv_pre_DEBIAS_",IDs[i],"_AUCs.csv",sep=""), row.names = FALSE)
write.csv(roc.df, paste("./output/",output,"/CSVs/builtenv_pre_DEBIAS_",IDs[i],"_ROCs.csv",sep=""), row.names = FALSE)
write.csv(feat_imp_df, paste("./output/",output,"/CSVs/builtenv_pre_DEBIAS_",IDs[i],"_featimp.csv",sep=""), row.names = FALSE)

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
pval.df <- data.frame(Study_ID = IDs[i], pval = for.pval, DEBIAS = FALSE)


if(file.exists(pval_filename)){
  all <- read.csv(pval_filename)
  all <- rbind(all, pval.df)
  write.csv(all, pval_filename, row.names = FALSE)
}else{(write.csv(pval.df, pval_filename, row.names = FALSE))}

print("done writing to CSVs")
print("done with all pre-DEBIAS")
