# If this ends up being a long code, maybe rename it to [James "Logan"] Howlett

# Data generation for the training of Cerebro.
# Burst relation based on the data infered from https://arxiv.org/pdf/2104.09536.pdf from figure 2 (p. 8)
# Most frequent cases for current SFR ranges from 1e-4 to 1e-2 Msol/yr/kpc^2 for 1e7 to 1e9 total mass of Msol/kpc^2

setwd("~/Documents/GitHub/LOGAN-SFH")
source("MUTANTS/LOGAN.R")
require(plyr)
library(ProSpect)
library(ggplot2)
library(plyr)         # splat(func)(c(var1, list_of_vars))
library(stringi)
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
'%!in%' <- function(x,y)!('%in%'(x,y))



# File should have the file extension '.rda'
outputFolder = "/Volumes/Elements/TrainingData/"
savefilename = "Argument_df.rda"
configFilename = "Data_Generation_Parameters.R"



# Load Parameters from the Config File.
source(configFilename)



generateDataFrameArguments <- function(Parameters,
                                       n.simul,
                                       speclib,
                                       save_path=NULL,
                                       progress_verbose=NULL,
                                       verbose=1){
  # Generates a DataFrame with all the parameters that are going to be used for the whole data.

  
  evaluateFunction <- function(x){
    # Evaluate a parameter that has a function
    additional_arguments = list(islog10=FALSE, min_val=NULL, max_val=NULL)
    for (n in names(additional_arguments)){
      if (n %in% names(x)){
        additional_arguments[n] = x[[n]]
        x <- x[-which(names(x)==n)]
      }
    }
    
    evaluated = splat(x$func)(c(1, x[-which(names(x)=="func")]))
    
    # Verify it is located within min/max
    if (!is.null(additional_arguments$min_val)){
      if (evaluated < additional_arguments$min_val){
        evaluated <- additional_arguments$min_val
      }
    }
    if (!is.null(additional_arguments$max_val)){
      if (evaluated > additional_arguments$max_val){
        evaluated <- additional_arguments$max_val
      }
    }
    
    # Evaluate exponent
    if (additional_arguments$islog10) {
      evaluated <- 10**evaluated
    }
    
    return(evaluated)
  }
  
  
  uniqueArguments <- function(Par, speclib){
    # This function returns a single set of parameters obtained from Parameters
    for (param in names(Par)){
      par_val = Par[[param]]
      
      # Functions
      if (is.list(par_val)){
        if ("func" %in% names(par_val)){
          # Evaluate the function
          tmp = evaluateFunction(par_val)
          Par[param] <- tmp
        } 
      # List/Sequence
      } else if (length(par_val) > 1){
        tmp = sample(par_val, 1)
        Par[param] <- tmp
      }
    }
  
    
    #### Adjust massfunction parameters ####
    # Get the arguments of the massfunction
    massfunc = Par[["massfunc"]]
    massfunc_args = Par[names(Par) %in% names(formals(massfunc))]
    agevec = speclib$Age
    ageweigths = speclib$AgeWeights
    ageweights_burst = ageweigths[agevec < Par[["mburstage"]]*1e9]
    
    # Calculate random probability of having burst
    if (!is.null(Par[["probburst"]])){
      if (Par[["probburst"]] >= 1){
        burst = TRUE
      } else if (0 <= Par[["probburst"]] && Par[["probburst"]] < 1){
        RndChanceBurst = runif(1, 0, 1)
        if (RndChanceBurst < Par[["probburst"]]){
          burst = TRUE
        } else {
          burst = FALSE
        }
      }
    } else {
      burst = TRUE
    }
    
    if (!burst){
      if (!is.null(Par[["mburst"]]) && is.null(Par[["burstfraction"]])){
        # Case A -> mburst is known, burstfraction is NULL
        Par[["mburst"]] = 0
      } else if (is.null(Par[["mburst"]]) && !is.null(Par[["burstfraction"]])){
        # Case A -> mburst is known, burstfraction is NULL
        Par[["burstfraction"]] = 0
      } else{
        cat("[PROBBURST] No specific case was found... This should be fixed!!!\n")
        print(Par)
        stop("NO CASE WAS FOUND!")
      }
    }
    
    if (!is.null(Par[["totalmass"]]) && !is.null(Par[["mburst"]]) && is.null(Par[["mSFR"]]) && is.null(Par[["burstfraction"]])){
      # Option A: totalmass and mburst known; mSFR=burstfraction=NULL
      if (verbose > 1) cat("Case A\n") # ToDo: Improve verbosity
      mass_total_noburst_SFR1 = sum(do.call("massfunc",
                                            c(list(agevec),
                                              massfunc_args[-which(names(massfunc_args)=="mSFR")],
                                              list(mSFR=1))
                                            ) * ageweigths)
      mass_total_burst = sum(Par[["mburst"]] * ageweights_burst)
      burstfraction = mass_total_burst / Par[["totalmass"]]
      if (burstfraction >= 1){
        if (verbose > 1) {
          cat(paste0("Case Abis: Mtotal_burst (according to given parameters) ",
                     "is greater than totalmass. Exclusively burst is going to be used.\n")) # ToDo: Improve verbosity
          }
        Par[["burstfraction"]] = 1
        Par[["mSFR"]] = 0
      } else {
        Par[["burstfraction"]] = burstfraction
        mSFR_scaled = (Par[["totalmass"]]*(1-burstfraction)) / mass_total_noburst_SFR1
        Par[["mSFR"]] = mSFR_scaled
      }
    } else if (!is.null(Par[["totalmass"]]) && is.null(Par[["mburst"]]) && is.null(Par[["mSFR"]]) && !is.null(Par[["burstfraction"]])){
      # Option B: totalmass and burstfraction known; mSFR=mburst=NULL
      if (verbose > 1) cat("Case B\n") # ToDo: Improve verbosity
      mass_total_noburst_SFR1 = sum(do.call("massfunc",
                                            c(list(agevec),
                                              massfunc_args[-which(names(massfunc_args)=="mSFR")],
                                              list(mSFR=1))
                                            ) * ageweigths)
      mass_total_burst = sum(Par[["burstfraction"]] * Par[["totalmass"]])
      mburst = mass_total_burst / sum(ageweights_burst)
      Par[["mburst"]] = mburst
      mSfR_scaled = (Par[["totalmass"]]*(1-Par[["burstfraction"]])) / mass_total_noburst_SFR1
      Par[["mSFR"]] = mSFR_scaled
    } else if (!is.null(Par[["totalmass"]]) && !is.null(Par[["mburst"]]) && !is.null(Par[["mSFR"]]) && is.null(Par[["burstfraction"]])){
      # Option C: totalmass, mSFR and mburst known; burstfraction=NULL
      if (verbose > 1) cat("Case C\n") # ToDo: Improve verbosity
      # Nothing to touch on
    } else if (is.null(Par[["totalmass"]]) && !is.null(Par[["mburst"]]) && !is.null(Par[["mSFR"]]) && is.null(Par[["burstfraction"]])){
      # Option D: mSFR and mburst known; totalmass=burstfraction=NULL
      if (verbose > 1) cat("Case D\n") # ToDo: Improve verbosity
      # Forcemass needs to be FALSE for SFHfunc
      Par[["totalmass"]] = FALSE
    } else {
      cat("[MASSFUNC] No specific case was found... This should be fixed!!!\n")
      print(Par)
      stop("NO CASE WAS FOUND!")
    }
    
    return(Par)
  }
  
  
  #### Initiate execution ####
  # If required parameters are not available, set the default ones.
  pnames = names(Parameters)
  def_parameters = list(RndSeed=sample(1:1e7, 1),
                        SNR=0,
                        probburst=runif(1, 0, 1),
                        totalmass=FALSE,
                        emission=TRUE,
                        emission_scale="SFR",
                        burstfraction=NULL,
                        filters="HST",
                        z=1e-4)
  for (name in names(def_parameters)){
    if (name %!in% pnames){
      Parameters[name] <- def_parameters[[name]]
    }
  }
  set.seed(Parameters$RndSeed)
  
  
  ##### Initialize data.frame with values that are going to differ with each case ####
  argument_names = c("mSFR")
  for (n in names(Parameters)){
    if (length(Parameters[[n]])> 1){
      argument_names <- c(argument_names, n)
    }
  }
  for (n in c("RndSeed", "probburst", "massfunc", "zfunc")){
    if (n %in% argument_names){
      argument_names <- argument_names[-which(argument_names == n)]
    }
  }
  
  df <- data.frame(matrix(0, nrow=n.simul, ncol=length(argument_names)))
  names(df) <- argument_names
  
  
  #### Iterate over n.simul and populate the  ####
  for (i in 1:n.simul){
    if (!is.null(progress_verbose)){
      if (i %% progress_verbose == 0){
        cat(paste0("Current iteration: ", i, "\n"))
      }
    }
    a <- uniqueArguments(Parameters, speclib)
    tmp <- data.frame(a[argument_names])
    df[i,] <- tmp
  }
  cat("FINISHED\n")
  print(object.size(df), units = "MB")
  if (!is.null(save_path)){
    cat(paste0("Saved file in ", save_path, " .\n"))
    save(df, file=save_path)
  }
  return(list(df=df, Parameters=Parameters))
}


drawSFHFromDataFrame <- function(df,
                                 func,
                                 speclib=NULL,
                                 agevec=NULL,
                                 specificCases=NULL,
                                 progress_verbose=1000,
                                 alpha_color=10,
                                 ...){
  # Draws the SFHs that are going to be used
  dots = list(...)
  
  if (is.null(speclib) && is.null(agevec)){
    stop("Either 'speclib' or 'agevec' needs to be given to the function")
  } else {
    if (is.null(agevec)){
      agevec <- speclib$Age
    }
  }
  
  if (!is.null(specificCases)){
    df <- df[specificCases, ]
  }
  
  # Evaluate the values of each of the SFH
  points_matrix <- matrix(0, nrow=dim(df)[1], ncol=length(agevec))
  args_func <- names(df) %in% names(formals(func))
  i <-  1
  cat("Painting...\n")
  while (i <= dim(df)[1]){
    points_matrix[i, ] <- do.call("func", c(list(age=agevec), df[i, args_func]))
    if (!is.null(progress_verbose)){
      if (i %% progress_verbose == 0){
        cat(paste0("Calculating: ", i, "\n"))
      }
    }
    i <- i + 1
  }
  
  
  # Populate dots with default variables in case they are not defined
  # This avoids multiple instances of the same argument, and the possibility to change them
  plot_def_arguments = list(type="l",
                            xlab="Age", ylab="SFR",
                            col=mycol <- rgb(0, 0, 255, max=255, alpha=alpha_color),
                            xlim=c(agevec[1],14e9),
                            ylim=c(0, max(points_matrix)),
                            log="x")
  for (def_arg in names(plot_def_arguments)){
    if (def_arg %!in% names(dots)){
      dots[def_arg] <- plot_def_arguments[def_arg]
    }
  }
  
  
  # Initialize plot
  do.call("plot", c(list(x=agevec, y=points_matrix[1,]), dots))
  
  # Iterate over each line and draw the data
  i <- 2
  col.lines <- dots[["col"]]
  while (i <= dim(points_matrix)[1]){
    do.call("lines", c(list(x=agevec, y=points_matrix[i, ], col=col.lines)))
    if (!is.null(progress_verbose)){
      if (i %% progress_verbose == 0){
        cat(paste0("Drawing: ", i, "\n"))
      }
    }
    i <- i + 1
  }
  return(points_matrix)
}


generateSpecFromDataFrame <- function(Parameters,
                                      df,
                                      verbose=1,
                                      stellpop="EMILESCombined",
                                      speclib=NULL,
                                      new_scale="defaultlog1",
                                      waveout=seq(4700, 9400, 1.25),
                                      ...){
  # Generate the Spectra from the matrix with the specific arguments
  
  #### Function definitioin ####
  insertNoise <- function(spec, SNRatio){
    # Each point in the spectrum will be modified
    # according to it's value and the SNRatio.
    spec$flux$flux = rnorm(length(spec$flux$flux), spec$flux$flux, spec$flux$flux/SNRatio)
    return(spec)
  }
  
  generateFilename <- function(rngLength=8){
    # Generates a filename
      now = Sys.time()
      filename = paste(gsub(" ", "T", gsub(":", "", gsub("-", "", now))),
                       stri_rand_strings(n=1, length=rngLength), sep="_")
  }
  
  #### Start of execution ####
  dots <- list(...)
  n.simul <- dim(df)[1]
  name.df <- names(df)
  par.names = names(Parameters)
  
  if (!is.null(speclib)){
    if (!is.null(stellpop) & stellpop != "EMILESCombined"){
      warning("stellpop and speclib do not match!")
    } else {
      # Return to one of the default accepted values for SFHfunc
      stellpop="EMILES"
    }
  } else {
    if (stellpop  == "EMILESCombined"){
      # Return to one of the default accepted values for SFHfunc
      stellpop = "EMILESCombined"
      # TODO: This will need to read global variable/env-variable or something similar.
      speclib = readRDS(file="EMILESData/EMILESCombined.rds")
    }
  }
  agevec = speclib$Age
  
  
  if (Parameters$filters == "HST"){
    filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
    filters <- list()
    for (filter in filtersHST) {
      # TODO: see speclib regarding location of data files.
      filters[[filter]] <-
        read.table(
          paste0("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat"),
          col.names = c("wave", "response")
        )
    }
    Parameters$filters <- filters
  }
  
  
  # Filter the dataframe according to where the info is going to go
  list.names.not.used = c()
  # massfunc
  massfunc.names = names(formals(Parameters$massfunc))
  massfunc_args = name.df[name.df %in% massfunc.names]
  list.names.not.used = c(list.names.not.used, massfunc.names[massfunc.names %!in% massfunc_args])
  # zfunc
  if (!is.function(Parameters$zfunc)){
    # This occurs if a value for Z was given.
    if (Parameters$zfunc > 1 || Parmeters$zfunc < 0){
      stop(paste0("Value fo Parameters$zfunc is not valid: ",
                  Parameters$zfunc, "\n"))
    }
  } else {
    zfunc.names = names(formals(Parameters$zfunc))
    zfunc_args = name.df[name.df %in% zfunc.names]
    list.names.not.used = c(list.names.not.used, zfunc.names[zfunc.names %!in% zfunc_args])
  }
  # SFHfunc
  SFHfunc.names = names(formals(SFHfunc))
  SFHfunc_args = name.df[name.df %in% SFHfunc.names]
  list.names.not.used = c(list.names.not.used, SFHfunc.names[SFHfunc.names %!in% SFHfunc_args])
  
  # Clean list of cases that are added manually
  list.values.considered <- c("massfunc", "...", "age", "forcemass", "Z",
                              "stellpop", "speclib", "filters", "emission",
                              "emission_scale")
  list.names.not.used = list.names.not.used[list.names.not.used %!in% list.values.considered]
  
  # Get constant values from Parameters that are also going to be inserted
  static.params.names = names(Parameters)[names(Parameters) %in% list.names.not.used]
  static.params.val = Parameters[static.params.names]
  
  
  
  #### Iterate over total number of cases ####
  # Initiate Matrix
  agevec_new = convertAgevecToOutputScale(agevec, agevec, new_scale = new_scale, return_scale = TRUE)
  numberColumnsIn = length(waveout) + length(Parameters$filters) + 1    # Add 1 for the ID
  numberColumnsLa = 2 * length(agevec_new$age) + 1           # Add 1 for the ID
  completeDataMatrixIn <- matrix(0, nrow=n.simul + 1, ncol=numberColumnsIn)
  completeDataMatrixLa <- matrix(0, nrow=n.simul + 1, ncol=numberColumnsLa)
  completeDataMatrixIn[1, ] = c(0, waveout, seq(1:length(filters)))
  completeDataMatrixLa[1, ] = c(0, agevec_new$age, agevec_new$age)
  
  
  i = 1
  while (i <= n.simul){
    
    # Evaluate elements that change each iteration
    # drop=FALSE keeps the name, even though there are none or just one column
    massfunc_args_val = df[i, massfunc_args, drop=FALSE]
    zfunc_args_val = df[i, zfunc_args, drop=FALSE]
    SFHfunc_args_val = df[i, SFHfunc_args, drop=FALSE]
    current_arguments = c(massfunc_args_val, zfunc_args_val,
                          SFHfunc_args_val, static.params.val)
    
    
    #### execute SHF ####
    spectraObject = do.call("SFHfunc", c(list(massfunc=Parameters$massfunc,
                                              forcemass=df[i, "totalmass"],
                                              Z=Parameters$zfunc,
                                              stellpop = stellpop,
                                              speclib = speclib,
                                              filters = Parameters[["filters"]],
                                              emission = Parameters$emission,
                                              emission_scale = Parameters$emission_scale
                                              ),
                                         current_arguments
                                         )
                            )
    a <- c(list(massfunc=Parameters$massfunc,
                forcemass=df[i, "totalmass"],
                Z=Parameters$zfunc,
                stellpop = stellpop,
                speclib = speclib,
                filters = Parameters[["filters"]],
                emission = Parameters$emission,
                emission_scale = Parameters$emission_scale),
           current_arguments)
    
    
    
    #### Postprocess Spectra ####
    # SFR
    tmp_output = convertAgevecToOutputScale(spectraObject$agevec, spectraObject$SFR)
    spectraObject$SFR = tmp_output$data
    # massvec
    tmp_output = convertAgevecToOutputScale(spectraObject$agevec, spectraObject$massvec)
    spectraObject$massvec = tmp_output$data
    # Zvec
    tmp_output = convertAgevecToOutputScale(spectraObject$agevec, spectraObject$Zvec)
    spectraObject$Zvec = tmp_output$data
    # agevec
    tmp_output = convertAgevecToOutputScale(spectraObject$agevec, spectraObject$agevec)
    spectraObject$agevec = suppressWarnings(tmp_output$data)    
    
    # Add Noise
    if (Parameters$SNR > 0){
      spectraObject = insertNoise(spectraObject, Parameters$SNR)
    }
    
    # Adjust Spectra to Wavelength (waveout)
    if (!is.null(waveout)){
      spectraObject$flux = interpolateToWaveout(
        lapply(spectraObject$flux["wave"], as.numeric)[[1]],
        lapply(spectraObject$flux["flux"], as.numeric)[[1]],
        waveout,
        returnList=TRUE)
    }
    
    # ID, spectra, Magnitudes, timer
    newRowIn <- c(i, spectraObject$flux$flux, spectraObject$out$out)
    newRowLa <- c(i, spectraObject$SFR, spectraObject$Zvec)
    
    # Add new row to the Matrix
    completeDataMatrixIn[i + 1, ] = newRowIn
    completeDataMatrixLa[i + 1, ] = newRowLa
    
    print(i)
    i = i + 1
  }
  
  #### Export File ####
  # Once calculations are over, export file
  outputfilename <- generateFilename()
  filterData = spectraObject$out
  UUIDs <-exportObjectsToSingleFITS(inputMatrix = completeDataMatrixIn,
                                    labelMatrix = completeDataMatrixLa,
                                    filename = outputfilename,
                                    foldername = "/Users/enrique/Documents/GitHub/LOGAN-SFH/KK",
                                    absolutePath = TRUE,
                                    filters = filterData,
                                    verbose=verbose
  )
  
  
  
  
  
}




#### Execute Code ####
output <- generateDataFrameArguments(Parameters=Parameters,
                                     n.simul=3,
                                     speclib=EMILESCombined,
                                     # save_path=file.path(outputFolder, savefilename),
                                     verbose=1,
                                     progress_verbose = 1000)
df = output[["df"]]
Parameters = output[["Parameters"]]

# point_matrix <- drawSFHFromDataFrame(df,
#                                      Parameters$massfunc,
#                                      agevec=EMILESCombined$Age,
#                                      log="xy",
#                                      ylim=c(1e-4, 15))

generateSpecFromDataFrame(Parameters, df)


# abline(v=EMILESCombined$Age)


