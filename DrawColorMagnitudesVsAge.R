# Draw Color magnitudes vs age (evolution of the magnitudes)

# Initialize variables
.pardefault <- par()
setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
library(RColorBrewer)
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

# Read EMILES Model and generate colorvector
zValues = EMILESCombined$Z
ageValues = EMILESCombined$Age
wave = EMILESCombined$Wave
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Dummy matrix with the magnitudes
data = matrix(NA, nrow=length(zValues), ncol=length(ageValues))

# List with the information
dataList = list()
for (filter in filtersHST){
  dataList[[filter]] = data
}


# Iterate over Z
for (zIdx in 1:length(zValues)){

  # Get corresponding spectra [age, wave]
  z.spect <- EMILESCombined$Zspec[zIdx][[1]]
  
  # Iterate over Ages
  for (ageIdx in 1:length(ageValues)){

    # Iterate over filters
    for (filter in filtersHST){
      
      # Evaluate bandpass for this filter, metallicity and age
      dataList[[filter]][zIdx, ageIdx] <- bandpass(wave, z.spect[ageIdx, ], filters[[filter]])
    }
  }
}



# Plot data
for (zIdx in 1:length(zValues)){
  
  pdf(file=sprintf("Plots/ColorvsAgevsZ/ColorEvolution_z%.4f.pdf", zValues[zIdx]))
  for (filterIdx in 1:length(filtersHST)){
    filter <- filtersHST[filterIdx]
    
    if (filterIdx == 1){
      plot(ageValues, dataList[[filter]][zIdx, ],
           type="l", col=col_vector[filterIdx], lwd=2,
           log="xy", xlab="Age", ylab="Flux", main=paste0("Flux Evolution over age for Z=", zValues[zIdx]))
    } else {
      lines(ageValues, dataList[[filter]][zIdx, ], lwd=2, col=col_vector[filterIdx])
    }
    
  }
  
  # legend
  legend(x = "topright",            # Position
         legend = filtersHST,       # Legend texts
         lty = 1,                   # Line types
         col = col_vector,          # Line colors
         lwd = 2)                   # Line width
  dev.off()
  
}


par(.pardefault)
