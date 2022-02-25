# Draw Color magnitudes vs age (evolution of the magnitudes)

# Initialize variables
.pardefault <- par()
setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
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

# Read EMILES Model
zValues = EMILESCombined$Z
ageValues = EMILESCombined$Age
wave = EMILESCombined$Wave

# Dummy matrix with the magnitudes
data = matrix(NA, nrow=length(zValues), ncol=length(ageValues))

# List with the information
dataList = list()
for (filter in filtersHST){
  dataList[[filter]] = data
}


# Iterate over Z
for (zIdx in 1:length(zValues)){
  z = zValues[zIdx]

  # Get corresponding spectra [age, wave]
  z.spect <- EMILESCombined$Zspec[zIdx][[1]]
  
  # Iterate over Ages
  for (ageIdx in 1:length(ageValues)){
    age = ageValues[ageIdx]

    # Iterate over filters
    for (filter in filtersHST){
      
      # Evaluate bandpass for this filter, metallicity and age
      dataList[[filter]][zIdx, ageIdx] <- bandpass(wave, z.spect[ageIdx, ], filters[[filter]])
    }
  }
}









par(.pardefault)
