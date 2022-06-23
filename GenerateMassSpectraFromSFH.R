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

# Define auxiliary parameters
waveout=seq(4700, 7500, 1.25)
filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filt_cenwave = c(2708.680, 3358.580, 4329.475, 5331.898, 8072.917)
filters <- list()
for (filter in filtersHST) {filters[[filter]] <- read.table(paste0("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat"), col.names = c("wave", "response"))}
rm(filtersHST, filter)


temp_folder = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder/"
files = list.files(temp_folder)
n_models = (length(files) - 9)/2
n_models = 1
cores = 4
start = 1
end = 8200




# input_spectra <- np$load(normalizePath(file.path(temp_folder, "input_spectra.npy")))
# input_magnitudes <- np$load(normalizePath(file.path(temp_folder, "input_magnitudes.npy")))
agevec <- np$load(normalizePath(file.path(temp_folder, "agevec.npy")))
# ageweight <- np$load(normalizePath(file.path(temp_folder, "ageweight.npy")))
# spectra_lambda <- np$load(normalizePath(file.path(temp_folder, "spectra_lambda.npy")))
# label_sfh <- np$load(normalizePath(file.path(temp_folder, "label_sfh.npy")))
# label_z <- np$load(normalizePath(file.path(temp_folder, "label_z.npy")))
# label_sfh_no_normalization <- np$load(normalizePath(file.path(temp_folder, "label_sfh_no_normalization.npy")))
# label_z_no_normalization <- np$load(normalizePath(file.path(temp_folder, "label_z_no_normalization.npy")))

predicted = list()
for (model_id_1based in 1:n_models){
    model_id_0based = model_id_1based - 1
    predicted[[model_id_0based + 1]] = list(sfh=np$load(normalizePath(file.path(temp_folder, paste0("predict_sfh_", model_id_0based, ".npy")))),
                                            z=np$load(normalizePath(file.path(temp_folder, paste0("predict_z_", model_id_0based, ".npy")))))
    # predict_sfh <- np$load(normalizePath(file.path(temp_folder, paste0("predict_sfh_", model_id_0based, ".npy"))))
    # predict_z <- np$load(normalizePath(file.path(temp_folder, paste0("predict_z_", model_id_0based, ".npy"))))
}

generate_spectra <- function(input.list){
    predicted_spectra = matrix(NA, nrow=n.points, ncol=length(waveout))
    predicted_magnitudes = matrix(NA, nrow=n.points, ncol=length(filters))
    for (id in input.list$start:input.list$stop){
        spectraObject = SFHfunc(
            massfunc = massfunc_custom,
            speclib = input.list$speclib,
            filters=input.list$filters,
            Z = Z_custom,
            magevec=input.list$agevec,
            msfh=input.list$predicted[[input.list$model]]$sfh[id, ],
            Zagevec=input.list$agevec,
            Zsfh=input.list$predicted[[input.list$model]]$z[id, ],
            emission=TRUE,
            emission_scale = "SFR"
        )
        spectraObject <- post_process_spectra(spectraObject, input.list$waveout)

        predicted_spectra[id,] <- spectraObject$flux$flux
        predicted_magnitudes[id, ] <- spectraObject$out$out

    }
    list(spectra=predicted_spectra, magnitudes=predicted_magnitudes)
}


# Prepare data for the paralellization
n.points = dim(predicted[[1]]$sfh)[1]
input_arguments <- list()
for (c in 1:cores){
    input_arguments[[c]] <- list(start=1 + (c-1) * end / cores,
                                 stop=c * end / cores,
                                 verbose=2,
                                 model=model_id_1based,
                                 agevec=agevec,
                                 waveout=waveout,
                                 predicted=predicted,
                                 filters=filters)
}

for (model_id_1based in 1:n_models){
    model_id_0based <- model_id_1based - 1
    predicted_spectra = matrix(NA, nrow=n.points, ncol=length(waveout))
    predicted_magnitudes = matrix(NA, nrow=n.points, ncol=length(filters))
    
    # Parallel
    start_time = Sys.time()
    cl <- makeCluster(2)
    clusterExport(cl, c("SFHfunc", "massfunc_custom", "Z_custom"), envir=environment())
    cat("Starting calculation of spectra...\n")
    output <- mclapply(X=input_arguments,
                       FUN=generate_spectra,
                       mc.cores=cores)
    stopCluster(cl)
    
    # Sort data into the correct placement
    predicted_spectra = matrix(NA, nrow=n.points, ncol=length(waveout))
    predicted_magnitudes = matrix(NA, nrow=n.points, ncol=length(filters))
    for (c in 1:cores){
        predicted_spectra[input_arguments[[c]]$start:input_arguments[[c]]$stop,] <- output[[c]]$spectra[input_arguments[[c]]$start:input_arguments[[c]]$stop,]
        predicted_magnitudes[input_arguments[[c]]$start:input_arguments[[c]]$stop,] <- output[[c]]$magnitudes[input_arguments[[c]]$start:input_arguments[[c]]$stop,]
    }
    
    end_time = Sys.time()
    print(end_time - start_time)
    
    spe_file <- suppressWarnings(normalizePath(file.path(temp_folder,
                                                        paste0("predict_spectra_",
                                                               model_id_0based,
                                                               ".npy"))))
    mag_file <- suppressWarnings(normalizePath(file.path(temp_folder,
                                                         paste0("predict_magnitudes_",
                                                                model_id_0based,
                                                                ".npy"))))
    cat("Saving spectra obtained with predicted sfh/z in: ",
        normalizePath(file.path(temp_folder)),
        "/predict_[spectra, magnitudes]_",
        model_id_0based,
        ".npy ...\n", sep="")
    np$save(spe_file, predicted_spectra)
    np$save(mag_file, predicted_magnitudes)
}


