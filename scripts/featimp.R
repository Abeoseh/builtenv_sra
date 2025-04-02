suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(biomformat))
#suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))


## example of running this in the terminal: Rscript ../scripts/featimp.R associated

setwd("C:/Users/brean/Downloads/masters/Fodor/builtenv_sra/output")
args <- commandArgs(trailingOnly = TRUE)
# input = "skin_floor"
input = args[1]
output = input

dir.create(file.path(input, "gini_v_gini"), showWarnings = FALSE)

print(paste("input/output folder:",input))

## list all the files for that comparison 
files = list.files(paste(input,"/CSVs",sep=""), "*featimp")

# lognorm <- read.csv(paste("./csv_files/",input,"/lognorm_data.csv",sep=""), check.names=FALSE)
phenos <- read.csv("../csv_files/phenotypes.csv") # for naming the graphs
IDs <- unlist(strsplit(files, "_"))[ c(4,9,14,19) ]
IDs = IDs[!is.na(IDs)]
print("IDs:")
print(IDs)

### Graph: ###
used_IDs = list()
for( i in 1:length(IDs) ){
  used_IDs = append(used_IDs, IDs[i])
  
  for( j in 2:length(IDs) ){
    if( IDs[i] != IDs[j] & !(IDs[j] %in% used_IDs)){
      print(paste(IDs[i], "vs", IDs[j]))
      df1 <- read.csv( paste(input,"/CSVs/" ,files[i], sep="") ) %>% select("feature", "MeanDecreaseGini")
      colnames(df1)[2] = "MeanDecreaseGini1"
      
      df2 <- read.csv( paste(input,"/CSVs/" ,files[j], sep="") ) %>% select("feature", "MeanDecreaseGini")
      colnames(df2)[2] = "MeanDecreaseGini2"
      
      df <- merge(df1, df2, by = "feature", all = TRUE)
      
      ## perform a correlation
      correlation = cor.test(df$MeanDecreaseGini1, df$MeanDecreaseGini2)
      pval = correlation$p.value 
      if (pval == 0){pval == 0} 
      else{pval = sprintf("%.2e", pval)}
      cor = correlation$estimate %>% signif(digits = 3)
      
      ## for naming the axes 
      phen1 <- filter(phenos, ID == IDs[i])
      phen1 <- phen1[1,2]
      
      phen2 <- filter(phenos, ID == IDs[j])
      phen2 <- phen2[1,2]
      

      png(paste("./",output,"/gini_v_gini/plot_",IDs[i],"v",IDs[j],".png",sep=""), width = 1050, height = 480)
      
      plot = ggplot(df, aes(x = MeanDecreaseGini1, y = MeanDecreaseGini2, label = feature)) +
        geom_point() +
        annotate("text", x=max(df$MeanDecreaseGini1), y = min(df$MeanDecreaseGini2)+.4, label = paste("cor:", cor)) +
        annotate("text", x=max(df$MeanDecreaseGini1), y = min(df$MeanDecreaseGini2), label = paste("p-value:", pval)) +
        labs(title = "Feature Importance vs Feature Importance plot", x = phen1, y = phen2) +
        theme(plot.title = element_text(size=22), axis.text=element_text(size=11),
              axis.title=element_text(size=15)) 
      
      
      print(plot)
      
      dev.off()
    }
  }
}

print("done")