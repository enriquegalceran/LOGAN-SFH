# Generate New Versions of EMILES 

setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)
library(FITSio)


source('Constants.R')
# source("FuncionesAuxiliaresYLibrerias.R")

######################################
# Vamos a ampliar ahora
metalicityMHtoZ = cbind(
  c(-2.32, -1.71, -1.31, -0.71, -0.40, 0.00, +0.22),
  c(0.0001, 0.0004, 0.001, 0.004, 0.008, 0.019, 0.03)
)
# YEMILESdirectory = "/Users/fta/Downloads/Young_EMILES/EMILES_SUPERYOUNG_KROUPA_UNIVERSAL"
YEMILESdirectory = "/Users/fta/Documents/RStudio/DataNewEMILES"
outmassesOld <-
  read.csv("out_mass_KU_PADOVA00.txt",
           sep = "",
           header = FALSE)
colnames(outmassesOld) <-
  c(
    "IMFType",
    "IMFSlope",
    "MH",
    "Age",
    "Mtotal",
    "MSrem",
    "MS",
    "Mrem",
    "MG",
    "MSremLv",
    "MSL",
    "Mv",
    "V13"
  )
outmassesYoung <-
  read.csv(file.path(YEMILESdirectory, "../out_mass_KU"),
           sep = "",
           header = FALSE)

colnames(outmassesYoung) <-
  c(
    "IMFType",
    "IMFSlope",
    "MH",
    "Age",
    "Mtotal",
    "MSrem",
    "MS",
    "Mrem",
    "MG",
    "MSremLv",
    "MSL",
    "Mv",
    "V13"
  )
outmassesCombined <-
  read.csv(file.path(YEMILESdirectory, "../outCombined"),
           sep = "",
           header = FALSE)
colnames(outmassesCombined) <-
  c(
    "IMFType",
    "IMFSlope",
    "MH",
    "Age",
    "Mtotal",
    "MSrem",
    "MS",
    "Mrem",
    "MG",
    "MSremLv",
    "MSL",
    "Mv",
    "V13"
  )

Ages = sort(unique(outmassesCombined$Age))
Z = sort(unique(outmassesCombined$MH))
MH = Z
for (i in 1:length(Z)) {
  Z[i]  = metalicityMHtoZ[which.min(abs(metalicityMHtoZ[, 1] - Z[i])), 2]
  MH[i] = metalicityMHtoZ[which.min(abs(metalicityMHtoZ[, 1] - MH[i])), 1]
}
Z = sort(unique(Z))
MH = sort(unique(MH))

# Rewrite MH to closest values for the filter
for (i in 1:length(outmassesCombined$MH)) {
  outmassesCombined$MH[i] = metalicityMHtoZ[which.min(abs(metalicityMHtoZ[, 1] -
                                                            outmassesCombined$MH[i])),
                                            1]
}

cat("Generating Zevo section\n")
Zevo = list()
for (metal in MH) {
  tmpOutZ = outmassesCombined[outmassesCombined$MH == metal,]
  tmpOutZ[, 13] = c(1, rep(0, length(tmpOutZ[, 1]) - 1))
  tmpdf = data.frame(tmpOutZ[, c(7, 9, 5, 13, 8)])
  colnames(tmpdf) <- c("SMstar", "SMgas", "SMtot", "SFR", "SMrem")
  Zevo = c(Zevo, list(tmpdf))
}


########################################
# Open directory and generate the list of files
files <- list.files(YEMILESdirectory)
numberFiles <- length(files)
lenNameFiles = nchar(files[1])
ZSol = 0.02

listOfMHs = c()
listOfAges = c()
for (file in files) {
  # Metallicity
  metal = substring(file, 9, 13)
  if (substring(metal, 1, 1) == "m")
    substring(metal, 1, 1) <- "-"
  else
    substring(metal, 1, 1) <- "+"
  if (metal == "+0.41")
    metal = "+0.22"
  listOfMHs = c(listOfMHs, as.numeric(metal))
  # Age
  age = substring(file, 15, 21)
  listOfAges = c(listOfAges, as.numeric(age) * 1e9)
}

# Generate list with Zs
metalsUnique = sort(unique(listOfMHs))
metalsUniqueZ = metalsUnique
for (i in 1:length(metalsUnique)) {
  metalsUniqueZ[i] = metalicityMHtoZ[which.min(abs(metalicityMHtoZ[, 1] - metalsUnique[i])), 2]
  metalsUnique[i] = metalicityMHtoZ[which.min(abs(metalicityMHtoZ[, 1] - metalsUnique[i])), 1]
}
metalsUnique = unique(metalsUnique)
metalsUniqueZ = unique(metalsUniqueZ)

listOfZs = listOfMHs
for (i in 1:length(listOfZs)) {
  listOfZs[i] = metalicityMHtoZ[which.min(abs(metalicityMHtoZ[, 1] - listOfZs[i])), 2]
  listOfMHs[i] = metalicityMHtoZ[which.min(abs(metalicityMHtoZ[, 1] - listOfMHs[i])), 1]
}

# List to identify if information is an ascii or a FITS
ascii = c()
for (name in files) {
  if (substring(name, nchar(name) - 4, nchar(name)) == ".fits") {
    ascii <- c(ascii, FALSE)
  } else {
    ascii <- c(ascii, TRUE)
  }
}

# Consolidate into a dataframe for easier filtering
tableFiles <-
  data.frame(files, listOfAges, listOfMHs, listOfZs, ascii)
tableFiles


# Generate Ages and AgesBin Vector.
# AgesBin[1] = 0, rest are values in between (linear).
# Last element of AgesBin = AgesBin[last] + diference between AgesBin[last-1] and AgesBin[last]
Ages = sort(unique(listOfAges))
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
endstr = "END                                                                             "
emptystr = "                                                                                "
Wave = NULL
# Sort by Metallicity
for (metal in metalsUnique) {
  # Matrix for this specific metallicity
  fprintf(
    "Generating data for metallicity: % 4.2f [M/H] - (%1d/%1d)\n",
    metal,
    match(metal, metalsUnique),
    length(metalsUnique),
    sep = ""
  )
  matTMP = matrix(0, nrow = length(Ages), ncol = 53689)
  mat_row = 1
  dimnames(matTMP)[[1]] <- c(rep("temp", length(Ages)))
  dimnames(matTMP)[[2]] <- NULL
  
  # Filter only relevant Metallicity
  tableTMP = tableFiles[tableFiles$listOfMHs == metal,]
  tableTMP = tableTMP[order(tableTMP$listOfAges),]    # Sort by ages
  for (file_i in 1:length(tableTMP$files)) {
    if (tableTMP$ascii[file_i] & tableTMP$listOfAges[file_i] == 63100000) {next}
    fprintf(
      "Generating matrix for: % 4.2f [M/H] - %5.2e Yrs - %62s - (%2d/%2d)\n",
      metal,
      tableTMP$listOfAges[file_i],
      tableTMP$files[file_i],
      file_i,
      length(tableTMP$files)
    )
    filepath = file.path(YEMILESdirectory, tableTMP$files[file_i])
    # Check if ascii file or FITS
    if (tableTMP$ascii[file_i]) {
      fileData <- read.csv(filepath)
      colnames(fileData) <- c("value")
      fileNElements = length(fileData$value)
      # Find first row with data
      searching = TRUE
      endposition = 1
      while (searching) {
        if (fileData$value[endposition] == endstr)
          searching = FALSE
        else
          endposition = endposition + 1
      }
      searching = TRUE
      startposition = endposition + 1
      while (searching) {
        if (fileData$value[startposition] == emptystr)
          startposition = startposition + 1
        else
          searching = FALSE
      }
      fileHeader = fileData$value[1:endposition]
      Data <- vector(mode = "numeric", length = fileNElements - startposition + 1)
      if (is.null(Wave)) {
        Wave <-
          vector(mode = "numeric", length = fileNElements - startposition + 1)
        for (i in 1:fileNElements - startposition + 1) {
          Wave[i] = as.numeric(substr(fileData$value[i + startposition - 1],  1, 19))
        }
      }
      for (i in 1:fileNElements - startposition + 1) {
        Data[i] = as.numeric(substr(fileData$value[i + startposition - 1], 20, 41))
      }
    } else {
      readData = readFITS(filepath)
      Data <- readData$imDat
    }
    
    # Add data to Matrix for specific metallicity
    matTMP[mat_row, ] <- Data
    mat_row = mat_row + 1
  }
  
  # Append matrix to list
  Zspec = c(Zspec, list(matTMP))
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
EMILESCombined <- list(
  Z = metalsUniqueZ,
  Age = Ages,
  AgeBins = AgesBin,
  AgeWeights = AgeWeights,
  Wave = Wave,
  Labels = Labels,
  Zspec = Zspec,
  Zevo = Zevo
)

################################
# Save EMILESCombined into a file
cat("Saving Stellar Population in", file.path(getwd(), "EMILESCombined.rds\n"))
saveRDS(EMILESCombined, "EMILESCombined.rds")

# Split File into multiple parts
for (i in 1:8){
  if (i == 7){
    for (j in 1:7){
      saveRDS(EMILESCombined[i][[1]][j], paste0("EMILES", i, "_", j, ".rds"))
    }
  } else {
    saveRDS(EMILESCombined[i], paste0("EMILES", i, ".rds"))
  }
}


##################################
# Verify both libraries are equal
# Rebuild
EMILESReconstr = list()
for (i in 1:6){
  EMILESReconstr = c(EMILESReconstr, readRDS(file=paste0("EMILES", i, ".rds")))
}
tmp = list()
for (j in 1:7){
  tmp = c(tmp, readRDS(file=paste0("EMILES7_", j, ".rds")))
}
EMILESReconstr = c(EMILESReconstr, list(Zspec=tmp))
EMILESReconstr = c(EMILESReconstr, readRDS(file="EMILES8.rds"))

# Verify Reconstructed is the same
for (i in 1:8){
  if (i == 7){
    for (j in 1:7){
      for (k in 1:60){
        if (!setequal(EMILESCombined[i][[1]][[j]][k,],
                      EMILESReconstr[i][[1]][[j]][k,])){
          cat("Error in: 7 - ", j, "-",  k, "\n")
        }
      }
    }
  } else {
    if (!setequal(EMILESCombined[i],
                  EMILESReconstr[i])){
      cat("Error in: ", i, "\n")
    }
  }
}
cat("Finished Comparing.\n")


# Load EMILESCombined from a file
# EMILESNew = readRDS(file="EMILESCombined.rds")

##################
# Release all temporary variables
# If you want to remove everything except final Emiles, uncomment next lines
# listOfVariables = ls()
# listOfVariables = c(listOfVariables[listOfVariables != "EMILESCombined"], "listOfVariables")
# rm(list = listOfVariables)

rm(
  a,
  age,
  Ages,
  AgesBin,
  AgeWeights,
  ascii,
  Data,
  EMILESReconstr,
  emptystr,
  endposition,
  endstr,
  file,
  file_i,
  fileData,
  fileHeader,
  fileNElements,
  filepath,
  files,
  i,
  Labels,
  lenNameFiles,
  listOfAges,
  listOfMHs,
  listOfZs,
  matTMP,
  mat_row,
  metal,
  metalicityMHtoZ,
  metalsUnique,
  metalsUniqueZ,
  MH,
  name,
  numberFiles,
  outmassesCombined,
  outmassesOld,
  outmassesYoung,
  readData,
  searching,
  startposition,
  tableFiles,
  tableTMP,
  tmpdf,
  tmpOutZ,
  Wave,
  YEMILESdirectory,
  Z,
  Zevo,
  ZSol,
  Zspec
)






