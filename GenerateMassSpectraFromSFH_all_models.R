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
HRpypop = readRDS(file="HrPyPopData/HRPyPop.rds")
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


temp_folder = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder_step2"
n_models = 7
cores = 7
start = 1
end = 1000
model_id_1based = 2



# input_spectra <- np$load(normalizePath(file.path(temp_folder, "input_spectra.npy")))
# input_magnitudes <- np$load(normalizePath(file.path(temp_folder, "input_magnitudes.npy")))
agevec <- np$load(normalizePath(file.path(temp_folder, "agevec.npy")))
# ageweight <- np$load(normalizePath(file.path(temp_folder, "ageweight.npy")))
# spectra_lambda <- np$load(normalizePath(file.path(temp_folder, "spectra_lambda.npy")))
# label_sfh <- np$load(normalizePath(file.path(temp_folder, "label_sfh.npy")))
# label_z <- np$load(normalizePath(file.path(temp_folder, "label_z.npy")))
# label_sfh_no_normalization <- np$load(normalizePath(file.path(temp_folder, "label_sfh_no_normalization.npy")))
# label_z_no_normalization <- np$load(normalizePath(file.path(temp_folder, "label_z_no_normalization.npy")))

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

generate_spectra <- function(input.list){
    predicted_spectra = matrix(NA, nrow=n.points, ncol=length(waveout))
    predicted_magnitudes = matrix(NA, nrow=n.points, ncol=length(filters))
    for (id in input.list$start:input.list$stop){
        spectraObject = SFHfunc(
            massfunc = massfunc_custom,
            speclib = EMILESCombined_reduced,
            filters=input.list$filters,
            Z = Z_custom,
            magevec=input.list$agevec,
            msfh=input.list$predicted$sfh[id, ],
            Zagevec=input.list$agevec,
            Zsfh=input.list$predicted$z[id, ],
            emission=TRUE,
            emission_scale = "SFR",
            z=1e-4,
            forcemass=input.list[["forcemass"]]
        )
        spectraObject <- post_process_spectra(spectraObject, input.list$waveout)
        
        predicted_spectra[id,] <- spectraObject$flux$flux
        predicted_magnitudes[id, ] <- spectraObject$out$out
        
    }
    list(spectra=predicted_spectra, magnitudes=predicted_magnitudes)
}

# Iterate over all models

for (model_id_1base in 0:n_models){
    model_id_0based = model_id_1base - 1
    cat(paste0("Iterating for model#", model_id_1base, "/", n_models, " ...\n"))
    
    if (model_id_1base == 0){
        predicted = list(sfh=np$load(normalizePath(file.path(temp_folder, "label_sfh_no_normalization.npy"))),
                         z=np$load(normalizePath(file.path(temp_folder, "label_z_no_normalization.npy"))))
    } else {
        predicted = list(sfh=np$load(normalizePath(file.path(temp_folder, paste0("predict_sfh_", model_id_0based, ".npy")))),
                         z=np$load(normalizePath(file.path(temp_folder, paste0("predict_z_",   model_id_0based, ".npy")))))
    }

    
    # Prepare data for the paralellization
    n.points = dim(predicted$sfh)[1]
    input_arguments <- list()
    for (c in 1:cores){
        input_arguments[[c]] <- list(start=ceiling(1 + (c-1) * end / cores),
                                     stop=min(c(ceiling(c * end / cores), end)),
                                     verbose=2,
                                     model=model_id_1based,
                                     agevec=agevec,
                                     waveout=waveout,
                                     predicted=predicted,
                                     filters=filters,
                                     forcemass=1e10)
    }
    for (i in 1:cores){
        cat("For core#", i, ": ", input_arguments[[i]]$start, "-", input_arguments[[i]]$stop, "\n", sep = "")
    }
    
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
    
    if (model_id_1base == 0){
        spe_file <- suppressWarnings(normalizePath(file.path(temp_folder,
                                                             "input_spectra_regenerated.npy")))
        mag_file <- suppressWarnings(normalizePath(file.path(temp_folder,
                                                             "input_magnitudes_regenerated.npy")))
    } else {
        spe_file <- suppressWarnings(normalizePath(file.path(temp_folder,
                                                             paste0("predict_spectra_",
                                                                    model_id_0based,
                                                                    ".npy"))))
        mag_file <- suppressWarnings(normalizePath(file.path(temp_folder,
                                                             paste0("predict_magnitudes_",
                                                                    model_id_0based,
                                                                    ".npy"))))
    }

    cat("Saving spectra obtained with predicted sfh/z in: ", spe_file, " ...\n", sep="")
    cat("Saving magnitudes obtained with predicted sfh/z in: ", mag_file, " ...\n", sep="")
    np$save(spe_file, predicted_spectra)
    np$save(mag_file, predicted_magnitudes)
}

