library(lubridate)
library(dplyr)

getwd()
setwd("C:/Users/brean/Downloads/masters/Fodor/builtenv_sra/")
getwd()



df <- read.csv("./PRJEB14474/PRJEB14474_SraRunTable_old.csv")
df$Date <- as.Date(df$Date, format="%m/%d/%Y")
df$month_year <- format(df$Date, "%Y/%m")

unique(df$month_year)

grouped <- group_by(df, month_year)  
View(count(grouped,sample_type))

df <- df[! ( df$month_year %in% c("0012/12", "0013/01") ),]

write.csv(df, "./PRJEB14474/PRJEB14474_SraRunTable.txt", row.names = F)
