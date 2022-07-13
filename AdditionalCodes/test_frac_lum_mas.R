#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

# Load parameters + libraries
original_parameters <- par()
library(ProSpect)
library(plotrix)
library(parallel)
setwd("~/Documents/GitHub/LOGAN-SFH")
source("MUTANTS/LOGAN-CLAWS.R")
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
EMILESCombined_reduced = readRDS(file="EMILESData/EMILESCombined_reduced.rds")
library(reticulate)
np <- import("numpy")
configFilename = "Data_Generation_Parameters.R"
source(configFilename)
temp_folder = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder/"


# Define auxiliary parameters
waveout=seq(4700, 7500, 1.25)
filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filt_cenwave = c(2708.680, 3358.580, 4329.475, 5331.898, 8072.917)
filters <- list()
for (filter in filtersHST) {filters[[filter]] <- read.table(paste0("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat"), col.names = c("wave", "response"))}
rm(filtersHST, filter)




model_id_1based = 3
model_id_0based = model_id_1based - 1




agevec <- np$load(normalizePath(file.path(temp_folder, "agevec.npy")))
predicted = list(sfh=np$load(normalizePath(file.path(temp_folder, paste0("predict_sfh_", model_id_0based, ".npy")))),
                 z=np$load(normalizePath(file.path(temp_folder, paste0("predict_z_",   model_id_0based, ".npy")))))


sfh = predicted$sfh[1, ]
Z = predicted$z[1,]

spectraObject = SFHfunc(
    massfunc = massfunc_custom,
    speclib = EMILESCombined,
    filters=filters,
    Z = Z_custom,
    magevec=agevec,
    msfh=sfh,
    Zagevec=agevec,
    Zsfh=Z,
    emission=TRUE,
    emission_scale = "SFR",
    z=1e-4
)
spectraObject = post_process_spectra(spectraObject, waveout)


plt.spectra <- function(spectraObject){
    
    plot(spectraObject$flux$wave, spectraObject$flux$flux, type="l", log="x")
    
}
plt.spectra(spectraObject)


# fraccmass = "Amount mass in that bin of total " [sum(fraccmass) == 1]

age = EMILESCombined$Age[1:(length(EMILESCombined$AgeWeights) - 2)]
ageweight = EMILESCombined$AgeWeights[1:(length(EMILESCombined$AgeWeights) - 2)]

fraccmass = ageweight/sum(ageweight)
plot(fraccmass)
print(sum(fraccmass))

# age_i * frac
age_frac = age * fraccmass
age_frac_sum = sum(age_frac)

# age_i * fracmas_i * (L/M)_i_z

age_l_m = list()
for (zi in 1:length(EMILESCombined$Z)){
    z_ = EMILESCombined$Z[[zi]]
    tmp = list(sum_lum_mas = rowSums(EMILESCombined$Zspec[[zi]], dims=1)[1:length(age_frac)])
    tmp = c(tmp, list(numerador = age_frac * tmp$sum_lum_mas))
    tmp = c(tmp, list(out=sum(tmp$numerador)/sum(tmp$sum_lum_mas)))
    age_l_m[[zi]] = tmp
}

names(age_l_m) <- EMILESCombined$Z

age_l_m_vector <- c()
for (i in 1:length(age_l_m)){
    age_l_m_vector = c(age_l_m_vector, age_l_m[[i]]$out)
}
plot(age_l_m_vector)

















