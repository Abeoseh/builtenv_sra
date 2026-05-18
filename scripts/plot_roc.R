library(dplyr)
getwd()
setwd("C:/Users/brean/Downloads/masters/Fodor/builtenv_sra")
getwd()

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
print("input: ")
print(args[1])

# input = "skin_floor"
# input = "skin_skinassoc_4"
# input = "associated_4"

input = "skin_skinassoc PRJEB11111/skin_skinassoc"

roc <- read.csv(paste("./output/",input,"/CSVs/builtenv_ROCs.csv",sep=""))
pval <- read.csv(paste("./output/",input,"/CSVs/builtenv_AUC_pvals.csv",sep=""))
auc <- read.csv(paste("./output/",input,"/CSVs/builtenv_AUCs.csv",sep=""))
study <- read.csv("./csv_files/phenotypes.csv")

IDs <- unique(roc$Study_ID)

for (ID in IDs){
  plot_title = study$Author[study$ID == ID]
  
  png(paste("./output/post_cluster/",input,"/pre_DEBIAS-M_RF_lognorm_ROC_",ID,".png",sep=""))
  
  ## set margins mar(bottom,left,top,right)
  # par("mar"=c(5,5,6,2))

  ## plot the abline
  plot(x = 1:0, y = 0:1, type = "l", axes=F, xlim=c(1,0), xlab = "Specificity", ylab = "Sensitivity")
  abline(a=1, b = -1 )
  ## set axes, box around graph,and title
  axis(2 ,at=seq(0,1, by=0.2), labels = paste(seq(0,1, by=.2)))
  axis(1 ,at=seq(1,0, by=-0.2), labels = paste(seq(1,0, by=-.2)))

  box(lwd = 1)
    
  title(paste("Training without: ", plot_title, sep=""), cex.main=1.5)

  
  ## plot permutaitions
  permutations = unique(roc$Permutation)
  for(permutation in 1:length(permutations)){
      current_roc <- filter(roc, Study_ID == !!ID & Permutation == permutation & DEBIAS == FALSE)
      # head(current_roc) %>% print()
      lines(current_roc$specificities, current_roc$sensitivities, col = alpha("black", 0.4) )
      # print("done")
    }
  
  ## plot actual data
  current_roc <- filter(roc, Study_ID == !!ID & Permutation == 0 & DEBIAS == FALSE)
  lines(current_roc$specificities, current_roc$sensitivities, col="red", lwd = 3)
  
  
  ## add permutations
  for.pval = pval$pval[pval$Study_ID == ID] %>% signif(3)
  for.auc = auc$AUC[auc$Study_ID == ID & auc$Permutation == FALSE] %>% round(3)
  ## add p-value to ROC plot
  if (for.pval == 0){text(x = 0.2, y = 0.01, label = "p<0.001", font = 2, col = "red", lwd = 3)
  }else{text(x = 0.2, y = 0.01, label = paste("p=", for.pval, sep=""), font = 2, col = "red", lwd = 3)}
  
  text(x = 0.2, y = 0.09, label = paste ("AUC:", for.auc), font = 2, col = "red", lwd = 3)
  
  dev.off()
  print(paste(ID,"done"))
}
