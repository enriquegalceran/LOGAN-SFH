# Draw Color magnitudes vs age (evolution of the magnitudes)


.pardefault <- par()
setwd("~/Documents/GitHub/LOGAN-SFH")
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")

filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filters <- list()
for (filter in filtersHST) {
  filters[[filter]] <-
    read.table(
      paste("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat", sep = ""),
      col.names = c("wave", "response")
    )
}


zValues = EMILESCombined$Z
age = EMILESCombined$Age

for (z_idx in 1:length(zValues)){
  
  
  
  
  
  
}









par(.pardefault)
