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

outputFolder = "/Volumes/Elements/TrainingData/"
configFilename = "Data_Generation_Parameters.R"

# Load Parameters from the Config File.
source(configFilename)

# set.seed(seed)

uniqueArguments <- function(Parameters){
  
  evaluateFunction <- function(x){
    if ("islog10" %in% names(x)){
      islog10 = x[["islog10"]]
      x <- x[-which(names(x)=="islog10")]
    } else {
      islog10 = FALSE
    }
    
    evaluated = splat(x$func)(c(1, x[-which(names(x)=="func")]))
    if (islog10) {
      evaluated <- 10**evaluated
    }
    return(evaluated)
  }
  
  # This function returns a single set of parameters to train give to SFHfunc
  
  # If required parameters are not available, set the default ones.
  pnames = names(Parameters)
  def_parameters = list(RndSeed=sample(1:1e7, 1),
                        probburst=runif(1, 0, 1),
                        totalmass=NULL,
                        div_massinburst_by_totalmass=FALSE,
                        z=1e-4)
  for (name in names(def_parameters)){
    if (name %!in% pnames){
      Parameters[name] <- def_parameters[[name]]
    }
  }
  
  
  # Evaluate all the parameters that are lists with "func" in one of the names.
  for (param in names(Parameters)){
    par_val = Parameters[[param]]
    if (is.list(par_val)){
      if ("func" %in% names(par_val)){
        # Evaluate the function
        tmp = evaluateFunction(par_val)
        Parameters[param] <- tmp
      } 
    } else if (length(par_val) > 1){
      tmp = sample(par_val, 1)
      Parameters[param] <- tmp
    }
  }
  
  
  
  
  
  return(Parameters)
}



a <- uniqueArguments(Parameters)

