"
  REMINDER FOR THE FUTURE:
  
  IF DEBUGGING IS GIVING OUT ERRORS (i.e. going to compiler::tryCmpFun), the 
  issue is because it is trying to debug the compiler (for a faster execution
  in the future). To fix this, the function needs to be 'precompiled':
    - set the working directory
    - Execute the source of LOGAN
    - Load libraries
    - COMMENT EVERY LIBRARY LINE OUT
    - Now it should work fine
  
"
setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
library(ggplot2)

EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")

outputFolder = "OutputsGenerated"
outputFolder = "/Volumes/Elements/Outputs"
outputFolder = "/Users/enrique/Documents/GitHub/LOGAN-SFH/OutputFolder"
absolutePath = TRUE
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
waveout = seq(4700, 9400, 1.25)
filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filters <- list()
for (filter in filtersHST) {
  filters[[filter]] <-
    read.table(
      paste("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat", sep = ""),
      col.names = c("wave", "response")
    )
}
  
# Stars = SFHfunc(
#   speclib = EMILESCombined,
#   filters = filters,
#   massfunc = massfunc_dtau,
#   mSFR = 10,
#   mpeak = 10,
#   emission = TRUE,
#   emission_scale = "SFR",
#   Z=Zfunc_massmap_box,
#   Zstart=1e-4,
#   Zfinal=0.1
# )

Stars = SFHfunc(
  speclib = EMILESCombined,
  filters = filters,
  massfunc = massfunc_snorm_burst,
  mSFR = 10,
  mpeak = 10,
  emission = TRUE,
  emission_scale = "SFR",
  mburst=5,
  Z=Zfunc_massmap_box,
  Zstart=1e-4,
  Zfinal=0.1
)
Stars2 = SFHfunc(
  speclib = EMILESCombined,
  filters = filters,
  massfunc = massfunc_snorm_burst,
  mSFR = 10,
  mpeak = 10,
  emission = TRUE,
  emission_scale = "SFR",
  mburst=0,
  Z=Zfunc_massmap_box,
  Zstart=1e-4,
  Zfinal=0.1
)



# plot(Stars)

