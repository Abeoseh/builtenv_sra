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
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(reprtree))

args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])
permutations = as.numeric(args[2])
input = args[3]
output = args[4]
batch = as.numeric(args[5])

lognorm_out <- read.csv(paste("./csv_files/",input,"/lognorm_data.csv",sep=""), colClasses = c("Phenotype" = "factor"))
lognorm_out <- lognorm_out[,c(1:3,6:length(lognorm_out))] # remove case and ID columns

phenos <- read.csv("./csv_files/phenotypes.csv") # for naming the graphs

## create AUC and ROC df
auc.df <- data.frame()

roc.df <- data.frame()

#### lognorm permutations testing as a loop ####

IDs <- unique(lognorm_out$Study_ID)
testing_ID = IDs[i]

random_forest_dataset <- filter(lognorm_out, Study_ID == testing_ID)
timepoints = unique(random_forest_dataset$timepoint)
print(dim(random_forest_dataset))
print("Study ID:")
print(testing_ID)
print(paste("timepoints in Study ID:", testing_ID))
print(timepoints)

if (length(timepoints) == 1){
  
  stop(paste("not enough timepoints in",testing_ID ,". The only timepoint is", timepoints))
}



# png(paste("./output/",output,"/plots/pre_DEBIAS-M_RF_lognorm_ROC_", IDs[i], ".png", sep=""))#, height = 24, width = 24)
# par(mar=c(3,3,1,0))#, mfrow=c(2,2))

delete_df <- data.frame()


for (training_timepoint in timepoints){
  training = filter(random_forest_dataset, timepoint == training_timepoint) %>% select(-c(sample_name, Study_ID,timepoint))
  print("__________________________________________________")
  print(paste("training ID:", training_timepoint))
  
  for (testing_timepoint in timepoints){
    if (testing_timepoint == training_timepoint & batch == 0){ 
      print("------------------------------------------------------------------------------------------------------------------")
      print(paste("testing ID:", testing_timepoint))
      testing = filter(random_forest_dataset, timepoint == testing_timepoint)  %>% select(-c(sample_name, Study_ID,timepoint))
      
      ## preform permutations
      training_permuted = training
      testing_permuted = testing
      for(j in 1:permutations){
        
        # set.seed(100) ## remove this line when removing if (batch == 0) 
        training_rows <- createDataPartition(y = training$Phenotype, p = .70, list = FALSE) ## remove this line when removing if (batch == 0) 
        testing_rows <- row.names(training)[! row.names(training) %in% training_rows] ## remove this line when removing if (batch == 0) 
        
        training_permuted <- training[training_rows,] ## remove this line when removing if (batch == 0) 
        testing_permuted <- testing[testing_rows, ] ## remove this line when removing if (batch == 0)... note training == testing  

        ## permutate training and testing labels
        training_permuted$Phenotype <- sample(training_permuted$Phenotype)
        testing_permuted$Phenotype <- sample(testing_permuted$Phenotype)
        ## do random forest
        RF_fit2 <- randomForest(Phenotype~., data = training_permuted)
        RF_pred <- predict(RF_fit2, testing_permuted, type = "prob")
        rf_roc <- roc(testing_permuted[,1], RF_pred[,1])
        
        # plot and add data to AUC and ROC dfs
        # if( j == 1){
        # p <- plot(rf_roc, print.auc=FALSE, add = FALSE)
        # }else{p <- plot(rf_roc, print.auc=FALSE, add = TRUE)}
        df = data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, AUC = auc(rf_roc), Permutation = TRUE, DEBIAS = "pre")
        auc.df <- rbind(auc.df, df)
        
        if (batch == 1){df <- data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, DEBIAS = "pre", Permutation = j, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)}else{
          df <- data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, DEBIAS = "pre", Permutation = (j+(permutations*(batch-1))), sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
        }
        roc.df <- rbind(roc.df, df)
        
        print(paste(j, " permutations done", sep=""))
      }
      ####
      #### Cross-fold Validation of data ####
      ####

      training_cv <- training
      cv_aucs <- c()

        ## if the batch is 0 and the training and testing dataset are the same do a 70-30 split with a 10-fold CV 
        for (cv in 1:10){ ## remove this line when removing "if (batch == 0)" 
          ## set seed and run rf
          # set.seed(100) ## remove this line when removing "if (batch == 0)" 
          
          training_rows <- createDataPartition(y = training_cv$Phenotype, p = .70, list = FALSE) ## remove this line when removing if (batch == 0) 
          testing_rows <- row.names(training_cv)[! row.names(training_cv) %in% training_rows] ## remove this line when removing if (batch == 0) 

          training <- training_cv[training_rows,] ## remove this line when removing "if (batch == 0)" 
          testing <- training_cv[testing_rows, ] ## remove this line when removing "if (batch == 0)" ... note training == testing  
          
          # set.seed(100)
          RF_fit <- randomForest(Phenotype~., data = training, importance=TRUE)
          
          # set.seed(100)    
          RF_pred <- predict(RF_fit, testing, type = "prob")
          
          ## compute ROC      
          rf_roc <- roc(testing[,1], RF_pred[,1])
          
          ## add AUC from ROC to AUC df
          # current_auc <- as.numeric(auc(rf_roc))
          # print(current_auc)
          cv_aucs <- append(as.numeric(auc(rf_roc)), cv_aucs)
          # print(cv_aucs)

          

          print("computed AUC of actual data")}
      
      

      df = data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, AUC = mean(as.numeric(cv_aucs)), Permutation = FALSE, DEBIAS = "pre") 
      auc.df <- rbind(auc.df, df)
      auc.df$AUC <- as.numeric(auc.df$AUC) # convert AUC in auc.df to numeric
      
      
      ## add a representative ROC to ROC df
      df <- data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, DEBIAS = "pre", Permutation = 0, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)        
      roc.df <- rbind(roc.df, df) 
      
      # print("-----------------model-----------------")
      # print(RF_fit)
      # 
      # delete_df1 <- data.frame(training_timepoint = rep(training_timepoint, length(dim(testing)[2])), testing_timepoint = rep(testing_timepoint, length(dim(testing)[2])), 
      #                          predictions = predict(RF_fit, testing), actual = testing$Phenotype)
      # delete_df <- rbind(delete_df, delete_df1)
      # 
      # png(paste("./output/",output,"/plots/decision_tree_training ",str_replace(training_timepoint, "/", "-"), " testing ", str_replace(training_timepoint, "/", "-"),".png", sep="")) 
      # reprtree:::plot.getTree(RF_fit)
      # dev.off()

    }else if (testing_timepoint != training_timepoint & batch != 0){ ## if testing and training timepoint are different
      print("------------------------------------------------------------------------------------------------------------------")
      print(paste("testing ID:", testing_timepoint))
      testing = filter(random_forest_dataset, timepoint == testing_timepoint)  %>% select(-c(sample_name, Study_ID,timepoint))
      
      ## preform permutations
      training_permuted = training
      testing_permuted = testing
      for(j in 1:permutations){
        
        set.seed(100) 
        ## permutate training and testing labels
        training_permuted$Phenotype <- sample(training_permuted$Phenotype)
        testing_permuted$Phenotype <- sample(testing_permuted$Phenotype)
        ## do random forest
        RF_fit2 <- randomForest(Phenotype~., data = training_permuted)
        RF_pred <- predict(RF_fit2, testing_permuted, type = "prob")
        rf_roc <- roc(testing_permuted[,1], RF_pred[,1])
        
        # plot and add data to AUC and ROC dfs
        # if( j == 1){
        # p <- plot(rf_roc, print.auc=FALSE, add = FALSE)
        # }else{p <- plot(rf_roc, print.auc=FALSE, add = TRUE)}
        df = data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, AUC = auc(rf_roc), Permutation = TRUE, DEBIAS = "pre")
        auc.df <- rbind(auc.df, df)
        
        if (batch == 1){df <- data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, DEBIAS = "pre", Permutation = j, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)}else{
          df <- data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, DEBIAS = "pre", Permutation = (j+(permutations*(batch-1))), sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)
        }
        roc.df <- rbind(roc.df, df)
        
        print(paste(j, " permutations done", sep=""))
      }
      if (batch == 1){
        ####
        #### predictions of actual data ####
        ####
        
        set.seed(100)
        RF_fit <- randomForest(Phenotype~., data = training, importance=TRUE)
        
        set.seed(100)
        RF_pred <- predict(RF_fit, testing, type = "prob")
        
        ## compute ROC      
        rf_roc <- roc(testing[,1], RF_pred[,1])
        
        df = data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, AUC = auc(rf_roc), Permutation = FALSE, DEBIAS = "pre") 
        auc.df <- rbind(auc.df, df)
        auc.df$AUC <- as.numeric(auc.df$AUC) # convert AUC in auc.df to numeric
        
        
        ## add ROC to ROC df
        df <- data.frame(Study_ID = IDs[i], training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, DEBIAS = "pre", Permutation = 0, sensitivities = rf_roc$sensitivities, specificities = rf_roc$specificities)        
        roc.df <- rbind(roc.df, df)
        
        # print("-----------------model-----------------")
        # print(RF_fit)
        # 
        # delete_df1 <- data.frame()
        # delete_df1 <- data.frame(training_timepoint = rep(training_timepoint, length(dim(testing)[2])), testing_timepoint = rep(testing_timepoint, length(dim(testing)[2])), 
        #                          
        #                          predictions = predict(RF_fit, testing), actual = testing$Phenotype)
        # delete_df1$actual_counts = NA
        # delete_df1$actual_counts[delete_df1$actual == 0] <- filter(count(training, Phenotype), Phenotype == 0)$n
        # delete_df1$actual_counts[delete_df1$actual == 1] <- filter(count(training, Phenotype), Phenotype == 1)$n
        # delete_df <- rbind(delete_df, delete_df1)
        # 
        # png(paste("./output/",output,"/plots/decision_tree_training ",str_replace(training_timepoint, "/", "-"), " testing ", str_replace(testing_timepoint, "/", "-"),".png", sep="")) 
        # reprtree:::plot.getTree(RF_fit)
        # dev.off()
        }}
    

    
}}



# if (batch == 0 | batch == 1) {
#   write.csv(delete_df, paste("./output/",output,"/CSVs/confusion_matrix_data_",batch,"_.csv",sep=""))
# 
# }




#### Write to CSV files ####
AUC_filename <- paste("./output/",output,"/CSVs/builtenv_AUCs.csv",sep="")
# pval_filename <- paste("./output/",output,"/CSVs/builtenv_AUC_pvals.csv",sep="")
ROC_filename <- paste("./output/",output,"/CSVs/builtenv_ROCs.csv",sep="")


write.csv(auc.df, paste("./output/",output,"/CSVs/builtenv_pre_DEBIAS_",IDs[i],"_batch_",batch,"_AUCs.csv",sep=""), row.names = FALSE)
write.csv(roc.df, paste("./output/",output,"/CSVs/builtenv_pre_DEBIAS_",IDs[i],"_batch_",batch,"_ROCs.csv",sep=""), row.names = FALSE)
# write.csv(feat_imp_df, paste("./output/",output,"/CSVs/builtenv_pre_DEBIAS_",IDs[i],"_featimp.csv",sep=""), row.names = FALSE)

## AUC
# if(file.exists(AUC_filename)){
#   all <- read.csv(AUC_filename)
#   all <- rbind(all, auc.df)
#   write.csv(all, AUC_filename, row.names = FALSE)
# }else{(write.csv(auc.df, AUC_filename, row.names = FALSE))}


## ROC
# if(file.exists(ROC_filename)){
#   all <- read.csv(ROC_filename)
#   all <- rbind(all, roc.df)
#   write.csv(all, ROC_filename, row.names = FALSE)
# }else{(write.csv(roc.df, ROC_filename, row.names = FALSE))}


## pvals
# pval.df <- data.frame(Study_ID = IDs[i], pval = for.pval, DEBIAS = FALSE)
# 
# 
# if(file.exists(pval_filename)){
#   all <- read.csv(pval_filename)
#   all <- rbind(all, pval.df)
#   write.csv(all, pval_filename, row.names = FALSE)
# }else{(write.csv(pval.df, pval_filename, row.names = FALSE))}

print("done writing to CSVs")
print("done with all pre-DEBIAS")
