# Generate New Versions of EMILES 

setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)
library(FITSio)
'%!in%' <- function(x,y)!('%in%'(x,y))

EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")

######################################
# Valores aproximados, de acuerdo con los otros valores (habria que haber hecho Z=10**MH)
metalicityMHtoZ = cbind(
  c(-1.74, -1.71, -1.33, -1.31, -0.72, -0.71,  -0.41, -0.40, 0.00, +0.22, +0.41),
  c(0.0004, 0.0004, 0.001, 0.001, 0.004, 0.004, 0.008, 0.008, 0.019, 0.03, 0.03)
)
# YEMILESdirectory = "/Users/fta/Downloads/Young_EMILES/EMILES_SUPERYOUNG_KROUPA_UNIVERSAL"
broadendirectory = "~/Documents/GitHub/LOGAN-SFH/miles_young_broaden"
files <- list.files(broadendirectory)
files <- files[-which(files == "readme.txt")]
files <- files[-which(files == "extra_data")]


EMILESRecortado = EMILESCombined
EMILESRecortado$Z = EMILESRecortado$Z[2:7]
EMILESRecortado$Zevo = EMILESRecortado$Zevo[2:7]
EMILESRecortado$Zspec = EMILESRecortado$Zspec[2:7]
EMILESRecortado$Age = EMILESRecortado$Age[1:58]

Z = EMILESRecortado$Z
Ages = EMILESRecortado$Age

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

EMILESRecortado$AgeBins = AgesBin
EMILESRecortado$AgeWeights = AgeWeights

# Fix Wave
new_wave = seq(from=2000, by=1.25, length.out=3400)
EMILESRecortado$Wave = new_wave


for (k in 1:length(EMILESRecortado$Zspec)){
  EMILESRecortado$Zspec[[k]] = matrix(NA, nrow=58, ncol=3400)
}

coverage = matrix(NA, nrow=58, ncol = 6)
for (file in files){
  
  # Read metallicty from filename
  metal = substring(file, 9, 13)
  if (substring(metal, 1, 1) == "m"){
    substring(metal, 1, 1) <- "-"
  } else {
    substring(metal, 1, 1) <- "+"
  }
  metal = as.numeric(metal)
  # Evaluate the metallicity from the dictionary and convert to z
  # (This step also fixes the issue with MH=+0.41 and small discrepancies
  # without needing to resort to min aligments)
  z = metalicityMHtoZ[which(metalicityMHtoZ[, 1] == metal), 2]
  
  # Read age from filename
  age = substring(file, 15, 21)
  age = as.numeric(age) * 1e9
  
  # Read data and select the matrix
  readData = readFITS(file.path(broadendirectory, file))
  Data <- readData$imDat
  # If array has two dimensions, squeeze to 1 dimension
  if (length(dim(Data)) == 2){
    Data = array(Data)
  }
  
  # evaluate which Z_idx and age_idx in our model
  Z_idx = which(Z == z)
  age_idx = which(Ages == age)
  
  # Substitute data in model
  EMILESRecortado$Zspec[[Z_idx]][age_idx, ] = Data[1:3400]
  coverage[age_idx, Z_idx] = 1
  
}

colnames(coverage) <- Z
rownames(coverage) <- Ages
no_NAs = sum(is.na(coverage))

obj_size <- object.size(EMILESCombined)
cat("Size of EMILESCombined: ", format(obj_size, "auto"), "\n")
obj_size <- object.size(EMILESRecortado)
cat("Size of EMILESRecortado: ", format(obj_size, "auto"), "\n")
cat("Total Amount of uncovered cases:", no_NAs, "\n")
cat("Saved EMILESRecortado in: ", normalizePath("EMILESRecortado.rds"), "\n")
saveRDS(EMILESRecortado, "EMILESRecortado.rds")

