library(ggplot2)
library(dplyr)
library(stringr)
getwd()
setwd("C:/Users/brean/Downloads/masters/Fodor/sink")
getwd()

input = "sink_nowater"

roc <- read.csv(paste("./output/",input,"/AUCs/builtenv_ROCs.csv",sep=""))
IDs <- unique(roc$Study_ID)
ID = IDs[1]


current_roc <- filter(roc, Study_ID == !!ID) %>% filter(Permutation == 0 & DEBIAS == FALSE)
png("./output/test.png")
# par("mar"=c(5,2,5,2))
par("mar"=c(5,2,8,2))
p <- plot(current_roc$specificities, current_roc$sensitivities, type="l", axes=F, xlim=c(1,0))
axis(2 ,at=seq(0,1, by=0.2), labels = paste(seq(0,1, by=.2)))
axis(1 ,at=seq(0,1, by=0.2), labels = paste(seq(1,0, by=-.2)))
title(str_wrap(paste("Training without: ", " European Kitchens: Moen et al. 2023", sep=""),30), line = + 8, cex.main=1.5)
print(p)
dev.off()







permutations = unique(roc$Permutation)
for(permutation in 1:length(permutations)){
    current_roc <- filter(roc, Study_ID == !!ID & Permutation == permutation & DEBIAS == FALSE)
    head(current_roc) %>% print()
    roc_plot <- roc_plot + 
      geom_line(current_roc, mapping = aes(x = specificities, y = sensitivities))
    print("done")
  }

print(roc_plot)
