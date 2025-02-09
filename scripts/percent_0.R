library(dplyr)

df = data.frame(x = c(0.5562311, 3.6971959, 3.2208731, 3.2208731, 1.0193071),
                y = c(NA, 0.1881033, 1.0193071, 0.1881033, 1.0193071),
                z = c(NA, NA, NA, 1.0193071, 0.1881033))


df = data.frame(x = c(0.5562311, 3.6971959, 3.2208731, 3.2208731, 1.0193071),
                y = c(0, 0.1881033, 1.0193071, 0.1881033, 1.0193071),
                z = c(0, NA, 0, 1.0193071, 0.1881033))
nrow(df) 

colSums(is.na(df))
colSums(df == 0)

percent = colSums(df == 0, na.rm = T)/nrow(df)
percent = percent[percent < 0.25] # gives the columns with an amount of NAs less than 25%



df %>% select(names(percent))


  
  
used_IDs = list()
for( ID1 in IDs ){
  used_IDs = append(used_IDs, ID1)
  for( ID2 in IDs ){
    if( ID1 != ID2 &  !(ID %in% used_IDs)){
      print(paste(ID1, "vs", ID2))}
	}
}
