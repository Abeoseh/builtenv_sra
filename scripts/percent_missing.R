#!/usr/bin/env/Rscript

.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

lognorm <- read.csv("csv_files/associated_2192/lognorm_data.csv", check.names = FALSE)
lognorm <- lognorm[6:length(lognorm)]
#### filtering rare taxa #### 

print("dim lognorm")
dim(lognorm)

# sum the amount of NA or 0 missing and divide by the amount of observations in that column (same number of observations in every column)
percent = colSums(lognorm == 0 | is.na(lognorm), na.rm = T)/nrow(lognorm) # percent is a named number
# make into a df
percent = data.frame(bacteria = names(percent), percent_missing = percent)

print("dim percent (number of rows in percent should be equal to number of columns in lognorm)")
dim(percent)

# percent = percent[percent < 0.25] # gives the columns with an amount of NAs less than 25%

write.csv(as.data.frame(percent), "./csv_files/test/percent_0s.csv")