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
                        probburst=runif(1, 0, 1),
                        totalmass=NULL,
                        burstfraction=NULL,
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
  return(df)
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


# drawSFHFromDataFrame(df, massfunc_snorm_burst, agevec=EMILESCombined$Age, main="Test1234")



generateSpecFromDataFrame <- function(Parameters,
                                      df){
  # Generate the Spectra from the matrix with the specific arguments
  
  
  
  
}




#### Execute Code ####
df <- generateDataFrameArguments(Parameters=Parameters,
                                 n.simul=1000,
                                 speclib=EMILESCombined,
                                 # save_path=file.path(outputFolder, savefilename),
                                 verbose=1,
                                 progress_verbose = 1000)

point_matrix <- drawSFHFromDataFrame(df,
                                     Parameters$massfunc,
                                     agevec=EMILESCombined$Age,
                                     log="xy",
                                     ylim=c(1e-4, 15))

# abline(v=EMILESCombined$Age)


