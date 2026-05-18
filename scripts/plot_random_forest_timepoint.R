.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
input = args[1]


#### read in csvs ####
# auc.df <- read.csv(paste("./output/",input,"/CSVs/","builtenv_ROCs.csv", IDs[i], ".png", sep=""))
# roc <- read.csv(paste("./output/",input,"/CSVs/","builtenv_ROCs.csv",sep=""))
study <- read.csv("./csv_files/phenotypes.csv")
amount_of_samples <- read.csv(paste("./output/", input,"/counts_after_filtering.csv", sep = ""))


#### make auc.df ####
auc_files = files = list.files(paste("./output/", input,"/CSVs",sep=""), "*_AUCs", full.names = TRUE)
auc_files[ grep("builtenv_AUCs", auc_files) ] <- NA
auc_files <- auc_files[!is.na(auc_files)]

# print(auc_files)

auc.df <- data.frame()
for (auc_file in auc_files){
  # print(auc_file)
  current_auc_df <- read.csv(auc_file)
  auc.df <- rbind(auc.df, current_auc_df)
}
write.csv(auc.df, paste("./output/", input,"/CSVs/builtenv_AUCs.csv",sep=""), row.names = F)




#### make roc df ####
roc_files = files = list.files(paste("./output/", input,"/CSVs",sep=""), "*_ROCs", full.names = TRUE)
roc_files[ grep("builtenv_ROCs", roc_files) ] <- NA
roc_files <- roc_files[!is.na(roc_files)]


roc <- data.frame()

for (roc_file in roc_files){
  # print(auc_file)
  current_roc_df <- read.csv(roc_file)
  roc <- rbind(roc, current_roc_df)
}
write.csv(roc, paste("./output/", input,"/CSVs/builtenv_ROCs.csv",sep=""), row.names = F)

IDs <- unique(roc$Study_ID)
print(paste("IDs:", IDs))


#### make p-val.df ####
pval.df <- data.frame()
for (ID in IDs){
  current_auc <- filter(auc.df, Study_ID == !!ID)

  training_timepoints = unique(auc.df$training_timepoint)
  testing_timepoints = unique(auc.df$testing_timepoint)
  
  # print(timepoints)
  for (training_timepoint in training_timepoints){
    for (testing_timepoint in testing_timepoints){
      if ( training_timepoint != testing_timepoint ){
        # print(timepoint)
        current_auc_df <- filter(current_auc, training_timepoint == !!training_timepoint & testing_timepoint == !!testing_timepoint)
        # print(current_auc_df)
    
        a <- current_auc_df[current_auc_df$Permutation == FALSE,]$AUC
    
        samp <- current_auc_df[current_auc_df$Permutation == TRUE,]$AUC
        for.pval = length(samp[samp >= a])/1000
        current_pval_df <- data.frame(Study_ID = ID, training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, pval = for.pval)
        pval.df <- rbind(pval.df, current_pval_df)
    
    	}
	}
  }
  print(paste("done with", ID))

}
print("done remaking AUC and ROC df")
write.csv(pval.df, paste("./output/",input,"/CSVs/builtenv_AUC_pvals.csv",sep=""), row.names=F)


#### make ROC curves ####
print("___________________________________________________")
print("plotting ROC curves")
for (ID in IDs){
  
  current_roc_studyid <- filter(roc, Study_ID == !!ID)
  
  training_timepoints = unique(current_roc_studyid$training_timepoint)
  testing_timepoints = unique(current_roc_studyid$testing_timepoint)
  print(training_timepoints)
  print(testing_timepoints)
  
    for (training_timepoint in training_timepoints){
      
      for (testing_timepoint in testing_timepoints){
      if (training_timepoint != testing_timepoint){
        current_roc = filter(current_roc_studyid, training_timepoint == !!training_timepoint & testing_timepoint ==  !!testing_timepoint)
        print(paste( "training:", unique(current_roc$training_timepoint) ))
        print(paste( "testing:", unique(current_roc$testing_timepoint) ))
        print("---------------------------------")
        
        training_timepoint_filename = gsub("/", "-", training_timepoint)
        testing_timepoint_filename = gsub("/", "-", testing_timepoint)
        
        png(paste("./output/",input,"/plots/training_",training_timepoint_filename,"_testing_", testing_timepoint_filename,"_pre_DEBIAS-M_RF_lognorm_ROC_", ID, ".png", sep=""))#, height = 24, width = 24)
        
        
      
        
        ## set margins mar(bottom,left,top,right)
        # par("mar"=c(5,5,6,2))
        
        ## plot the abline
        plot(x = 1:0, y = 0:1, type = "l", axes=F, xlim=c(1,0), xlab = "Specificity", ylab = "Sensitivity")
        abline(a=1, b = -1 )
        ## set axes, box around graph,and title
        axis(2 ,at=seq(0,1, by=0.2), labels = paste(seq(0,1, by=.2)))
        axis(1 ,at=seq(1,0, by=-0.2), labels = paste(seq(1,0, by=-.2)))
        
        box(lwd = 1)
        plot_title = study$Author[study$ID == ID]
    		
    			
    
    		## get number of training samples
        training_n = filter(amount_of_samples, timepoint == training_timepoint)$n	
        
        ## get number of testing samples
        testing_n = filter(amount_of_samples, timepoint == testing_timepoint)$n	
        
        title(str_wrap( paste("Train on: ", training_timepoint," (",training_n[2], ",", training_n[1],") Test on: ", testing_timepoint," (",testing_n[2], ",", testing_n[1],") for ", plot_title, sep=""), 40), cex.main=1.5)
        
        
        ## plot permutaitions
        # print("done 1")
        permutations = unique(roc$Permutation)
        # print("done 2")
        # print(permutations)
        for(permutation in 1:length(permutations)){
          permutation_roc <- filter(current_roc, Permutation == permutation)
          # print("done 3")
          # head(current_roc) %>% print()
          lines(permutation_roc$specificities, permutation_roc$sensitivities, col = alpha("black", 0.4) )
          # print("done")
        }
        
        ## plot actual data
        actual_data_roc <- filter(current_roc, Permutation == 0)
        # print("done 4")
        lines(actual_data_roc$specificities, actual_data_roc$sensitivities, col="red", lwd = 3)
        
        
        ## add permutations
        for.pval = pval.df$pval[pval.df$Study_ID == ID & pval.df$training_timepoint == training_timepoint & pval.df$testing_timepoint == testing_timepoint] %>% signif(3)
        # print(for.pval)
        # print("done 5")
        for.auc = auc.df$AUC[auc.df$Study_ID == ID & auc.df$Permutation == FALSE & auc.df$training_timepoint == training_timepoint & auc.df$testing_timepoint == testing_timepoint] %>% round(3)
        # print("done 6")
        ## add p-value to ROC plot
        if (for.pval == 0){text(x = 0.2, y = 0.01, label = "p<0.001", font = 2, col = "red", lwd = 3)
        }else{text(x = 0.2, y = 0.01, label = paste("p=", for.pval, sep=""), font = 2, col = "red", lwd = 3)}
        
        text(x = 0.2, y = 0.09, label = paste("AUC:", for.auc), font = 2, col = "red", lwd = 3)
        
        dev.off()
        
      
    
      }}}

  print(paste(ID,"done"))
  print("___________________________________________")
}

print("script complete")


