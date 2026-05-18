.libPaths( c( .libPaths(), "~/my_R_libs") )

suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(randomForest))
# suppressPackageStartupMessages(library(pROC))
# suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library(stringr))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(openxlsx))


#### getwd
getwd()
setwd("C:/Users/brean/Downloads/masters/Fodor/builtenv_sra")
getwd()

args <- commandArgs(trailingOnly = TRUE)
input = args[1]

# input  = "timepoint/skin_floor_timepoint"
# input  = "timepoint/associated_timepoint"
# input  = "timepoint/skinassoc_timepoint"

#### read in AUCs ####
study_names <- read.csv("./csv_files/phenotypes.csv") 
auc_files = list.files(paste("./output/", input,"/CSVs",sep=""), "*_AUCs", full.names = TRUE)
auc_files[ grep("builtenv_AUCs", auc_files) ] <- NA
auc_files <- auc_files[!is.na(auc_files)]


auc.df <- data.frame()
for (auc_file in auc_files){
  # print(auc_file)
  current_auc_df <- read.csv(auc_file)
  auc.df <- rbind(auc.df, current_auc_df)
}
# write.csv(auc.df, paste("./output/", input,"/CSVs/builtenv_AUCs.csv",sep=""), row.names = F)


IDs <- unique(auc.df$Study_ID)
print(paste("IDs:", IDs))

pval.df <- data.frame()
for (ID in IDs){
  current_auc <- filter(auc.df, Study_ID == !!ID)
  
  training_timepoints = unique(auc.df$training_timepoint)
  testing_timepoints = unique(auc.df$testing_timepoint)
  
  # print(timepoints)
  for (training_timepoint in training_timepoints){
    for (testing_timepoint in testing_timepoints){
      #if ( training_timepoint != testing_timepoint ){
        # print(timepoint)
        current_auc_df <- filter(current_auc, training_timepoint == !!training_timepoint & testing_timepoint == !!testing_timepoint)
        # print(current_auc_df)
        
        a <- current_auc_df[current_auc_df$Permutation == FALSE,]$AUC
        
        samp <- current_auc_df[current_auc_df$Permutation == TRUE,]$AUC
        if (training_timepoint != testing_timepoint){for.pval = length(samp[samp >= a])/1000}else{for.pval = length(samp[samp >= a])/100}
        
        current_pval_df <- data.frame(Study_ID = ID, training_timepoint = training_timepoint, testing_timepoint = testing_timepoint, pval = for.pval)
        pval.df <- rbind(pval.df, current_pval_df)
        
        #}
      }
  }
  print(paste("done with", ID))
  
}


# auc.df %>% filter(Permutation == F)

#### Beginning for loop to make excels of AUCs ####
for (ID in IDs){
    
    
    
    auc_pval_df <- merge(filter(auc.df, Permutation == F), pval.df, by = c("Study_ID", "training_timepoint", "testing_timepoint")) %>% filter(Study_ID == !!ID) %>%
      merge(study_names, all.x = T, by.x = "Study_ID", by.y = "ID") %>% relocate(Author) %>% select(!Study_ID)
    
    
    
    
    #### make a df of aucs in matrix form #### 
    ## testing dates are column names while training dates are row names
    timepoints <- unique(auc_pval_df$training_timepoint)
    amount_of_timepoints <- length(timepoints)
                                   
    auc_matrix <- setNames(data.frame(matrix(ncol = amount_of_timepoints, nrow = amount_of_timepoints), row.names = timepoints), timepoints)
    auc_matrix
    
    for (column_timepoint in timepoints){
      for (row_timepoint in timepoints){
        auc_matrix[row_timepoint, column_timepoint] <- auc_pval_df$AUC[auc_pval_df$testing_timepoint == column_timepoint & auc_pval_df$training_timepoint == row_timepoint]
        
      }
      
    }
    
    #### color AUC based off significance #### 
    
    
    # 3. Set up workbook
    wb <- createWorkbook()
    addWorksheet(wb, "FormattedTable")
    
    # 4. Write data (add row names as first column)
    auc_matrix <- cbind(Training = rownames(auc_matrix), auc_matrix)
    writeData(wb, "FormattedTable", auc_matrix, startRow = 2)
    
    ## add merged testing label
    mergeCells(wb, "FormattedTable", cols = 1:amount_of_timepoints+1, rows = 1)
    writeData(wb, "FormattedTable", "Testing", startRow = 1, startCol = 2)
    writeData(wb, "FormattedTable", sapply(strsplit(unique(auc_pval_df$Author), ":" ), `[`, 1), startRow = 1, startCol = 1) ## add the study to the table
    
    addStyle(
      wb, "FormattedTable",
      style = createStyle(halign = "center", fontSize = 11),
      rows = 1, cols = 1:amount_of_timepoints+1,
      gridExpand = TRUE
    ) ## center "timepoints" in merged cells
    
    

    # 5. Create a style (e.g., fill for requirements > 10)
    pval_significant_but_not_0 <- createStyle(fgFill = "#a60000", fontColour = "#ff5b5b") 
    # pval_0 <- createStyle(fgFill = "lightgreen")
    # pval_significant_but_not_0_text <- createStyle(fontColour = "red")

    number_style <- createStyle(numFmt = "#,##0")
    
    # 6. Loop over auc_pval_df to apply formatting
    for (i in seq_len(nrow(auc_pval_df))) {
      row_lab <- as.character(auc_pval_df$training_timepoint[i])
      col_lab <- as.character(auc_pval_df$testing_timepoint[i])
      pval <- auc_pval_df$pval[i]
      
      # Find matching row/col in auc_matrix
      row_index <- which(auc_matrix$Training == row_lab) + 2  # +1 because of header
      col_index <- which(colnames(auc_matrix) == col_lab)
      
      if (pval < 0.05) {
        addStyle(
          wb, "FormattedTable", pval_significant_but_not_0,
          rows = row_index, cols = col_index, stack = TRUE
        )}
        # else if (pval == 0){
        #   addStyle(
        #     wb, "FormattedTable", pval_0,
        #     rows = row_index, cols = col_index, stack = TRUE
        #   )}
      

    }

    
    # 7. Save the workbook
    saveWorkbook(wb, paste("./output/", input,"/", ID, "_colored_AUC.xlsx", sep=""), overwrite = TRUE)

}








