# Generate New Versions of EMILES 

setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)
library(FITSio)
library(stringr)
'%!in%' <- function(x,y)!('%in%'(x,y))
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
# source("FuncionesAuxiliaresYLibrerias.R")


HRDirectory = "/Volumes/Elements/HRpypopstar/sp"
HRDirectory_fractions = "/Volumes/Elements/HRpypopstar/remaining_stellar_mass_fractions"
HRSaving_directory = "/Volumes/Elements/HRpypopstar"
HRSaving_directory_split = "/Volumes/Elements/HRpypopstar/split"
IMF_used = "CHA"
list_of_all_files = list.files(HRDirectory)
list_files_IMF = list_of_all_files[which(substr(list_of_all_files, 10, 12) == IMF_used)]
numberFiles <- length(list_files_IMF)


# Calculate unique values
list_possibleZ = sort(unique(as.numeric(substr(list_files_IMF, 15, 19))))
list_possibleZs = sort(unique(substr(list_files_IMF, 15, 19)))

list_possiblelogt = sort(unique(as.numeric(substr(list_files_IMF, 25, 29))))
list_possiblelogts = sort(unique(substr(list_files_IMF, 25, 29)))

list_possibleAge = 10^list_possiblelogt


# Table with all the files
tableFiles <-
  data.frame(list_files_IMF, as.numeric(substr(list_files_IMF, 15, 19)), as.numeric(substr(list_files_IMF, 25, 29)))
colnames(tableFiles) <- c("filename", "Z", "logt")


# Read first element
data_1 <- read.table(file.path(HRDirectory,list_files_IMF[1]), sep = "", header = F, na.strings ="", stringsAsFactors= F)
length_spectra = dim(data_1)[1]
Wave = data_1[,1]



# Generate Ages and AgesBin Vector.
# AgesBin[1] = 0, rest are values in between (linear).
# Last element of AgesBin = AgesBin[last] + diference between AgesBin[last-1] and AgesBin[last]
Ages = sort(unique(list_possibleAge))
AgesBin = numeric(length(Ages) + 1)
for (a in 1:(length(Ages) - 1)) {
  AgesBin[a + 1] = (Ages[a] + Ages[a + 1]) / 2
}
AgesBin[length(AgesBin)] = Ages[length(Ages)] + AgesBin[length(AgesBin)] -
  AgesBin[length(AgesBin) - 1]
# AgeWeights is the total time in years for the stellar population bins
AgeWeights = numeric(length(Ages))
for (i in 1:length(AgeWeights)){
  AgeWeights[i] = AgesBin[i + 1] - AgesBin[i]
}




Zspec = list()
# Sort by Metallicity
for (metal in list_possibleZ) {
  # Matrix for this specific metallicity
  fprintf(
    "Generating data for metallicity: %4.3f Z - (%1d/%1d)\n",
    metal,
    match(metal, list_possibleZ),
    length(list_possibleZ)
  )
  matTMP = matrix(0, nrow = length(Ages), ncol = length_spectra)
  mat_row = 1
  dimnames(matTMP)[[1]] <- c(rep("temp", length(Ages)))
  dimnames(matTMP)[[2]] <- NULL
  
  # Filter only relevant Metallicity
  tableTMP = tableFiles[tableFiles$Z == metal,]
  tableTMP = tableTMP[order(tableTMP$logt),]    # Sort by ages
  for (file_i in 1:length(tableTMP$filename)) {
    fprintf(
      "Generating matrix for: % 4.2f Z - %5.2e logYrs - %s - (%2d/%2d)\n",
      metal,
      tableTMP$logt[file_i],
      tableTMP$filename[file_i],
      file_i,
      length(tableTMP$filename)
    )
    filepath = file.path(HRDirectory, tableTMP$filename[file_i])
    readData = read.table(filepath, sep = "", header = F, na.strings ="", stringsAsFactors= F)
    Data <- readData$V2
    
    # Add data to Matrix for specific metallicity
    matTMP[mat_row, ] <- Data
    mat_row = mat_row + 1
  }
  
  # Append matrix to list
  Zspec = c(Zspec, list(matTMP))
}





## ToDo: Zevo
filename = file.path(HRDirectory_fractions,
                     paste0("remaining_stellar_mass_fraction_", IMF_used, ".dat"))
outmasses = read.table(filename, sep = "", header = T, na.strings ="", stringsAsFactors= F)

Zevo = list()
for (metal in list_possibleZ) {
  tmpOutZ = outmasses[outmasses$Z == metal,]
  tmpOutZ[, 4] = c(1, rep(0, length(tmpOutZ[, 1]) - 1))
  tmpOutZ[, 5] = rep(1, length(tmpOutZ[,1]))
  tmpOutZ[, 6] = 1 - tmpOutZ[, 3]
  tmpOutZ[, 7] = rep(0, length(tmpOutZ[,1]))
  tmpdf = data.frame(tmpOutZ[, c(3, 6, 5, 4, 7)])
  colnames(tmpdf) <- c("SMstar", "SMgas", "SMtot", "SFR", "SMrem")
  Zevo = c(Zevo, list(tmpdf))
}


##################################################
# Paste everything together
# Last missing variable
Labels <- list(
  Zlab = "Metallicity",
  Agelab = "Time since ZAM / Yrs",
  Wavelab = "Wavelength / Ang",
  Lumlab = "Lsun / Ang (for 1 Msun SF)",
  LumAgeLab = "Lsun / Ang (for 1 Msun/Yr SFR)"
)

cat("Compiling final object\n")
HRPyPop <- list(
  Z = sort(unique(list_possibleZ)),
  Age = Ages,
  AgeBins = AgesBin,
  AgeWeights = AgeWeights,
  Wave = Wave,
  Labels = Labels,
  Zspec = Zspec,
  Zevo = Zevo
)

outputname <- paste0("HRPyPop_", IMF_used, ".rds")
cat("Saving HRPyPop in", file.path(HRSaving_directory, outputname),"\n")
saveRDS(HRPyPop, file.path(HRSaving_directory, outputname))


# Split File into multiple parts
cat("Saving HRPyPop split into multiple files in", file.path(HRSaving_directory, outputname),"\n")
for (i in 1:8){
  if (i == 7){
    for (j in 1:4){
      saveRDS(EMILESCombined[i][[1]][j], file.path(HRSaving_directory_split,
                                                   paste0("HRPyPop_", IMF_used, "_", i, "_", j, ".rds")))
    }
  } else {
    saveRDS(EMILESCombined[i], file.path(HRSaving_directory_split, paste0("HRPyPop_", IMF_used, "_", i, ".rds")))
  }
}



HRPypopReconstr= list()
for (i in 1:6){
  HRPypopReconstr = c(HRPypopReconstr, readRDS(file=file.path(HRSaving_directory_split, paste0("HRPyPop_", IMF_used, "_", i, ".rds"))))
}
tmp = list()
for (j in 1:4){
  tmp = c(tmp, readRDS(file=file.path(HRSaving_directory_split, paste0("HRPyPop_", IMF_used, "_7_", j, ".rds"))))
}
HRPypopReconstr = c(HRPypopReconstr, list(Zspec=tmp))
HRPypopReconstr = c(HRPypopReconstr, readRDS(file=file.path(HRSaving_directory_split, paste0("HRPyPop_", IMF_used, "_8.rds"))))












