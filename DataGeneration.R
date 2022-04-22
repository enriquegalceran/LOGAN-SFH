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

generateDataFrameArguments <- function(Parameters, n.simul, speclib){
  # Generates a DataFrame with all the parameters that are going to be used for the whole data.

  
  evaluateFunction <- function(x){
    # Evaluate a parameter that has a function
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
    ageweights_burst = ageweights[agevec < Par[["mburstage"]]*1e9]
    
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
    print(RndChanceBurst)
    print(Par$probburst)
    print(burst)
    if (!is.null(Par[["totalmass"]]) && !is.null(Par[["mburst"]]) && is.null(Par[["mSFR"]]) && is.null(Par[["burstfraction"]])){
      # Option A: totalmass and mburst known; mSFR=burstfraction=NULL
      cat("Case A\n") # ToDo: Improve verbosity
      mass_total_noburst_SFR1 = sum(do.call("massfunc", c(list(agevec), massfunc_args[-which(names(massfunc_args)=="mSFR")], list(mSFR=1))) * ageweigths)
      mass_total_burst = sum(Par[["mburst"]] * ageweights_burst)
      burstfraction = mass_total_burst / Par[["totalmass"]]
      if (burstfraction >= 1){
        cat("Case Abis: Mtotal_burst (according to given parameters) is greater than totalmass. Exclusively burst is going to be used.\n") # ToDo: Improve verbosity
        Par[["burstfraction"]] = 1
        Par[["mSFR"]] = 0
      } else {
        Par[["burstfraction"]] = burstfraction
        mSFR_scaled = (Par[["totalmass"]]*(1-burstfraction)) / mass_total_noburst_SFR1
        Par[["mSFR"]] = mSFR_scaled
      }
    } else if (!is.null(Par[["totalmass"]]) && is.null(Par[["mburst"]]) && is.null(Par[["mSFR"]]) && !is.null(Par[["burstfraction"]])){
      # Option B: totalmass and burstfraction known; mSFR=mburst=NULL
      cat("Case B\n") # ToDo: Improve verbosity
      mass_total_noburst_SFR1 = sum(do.call("massfunc", c(list(agevec), massfunc_args[-which(names(massfunc_args)=="mSFR")], list(mSFR=1))) * ageweigths)
      mass_total_burst = sum(Par[["burstfraction"]] * Par[["totalmass"]])
      mburst = mass_total_burst / sum(ageweights_burst)
      Par[["mburst"]] = mburst
      mSfR_scaled = (Par[["totalmass"]]*(1-Par[["burstfraction"]])) / mass_total_noburst_SFR1
      Par[["mSFR"]] = mSFR_scaled
    } else if (!is.null(Par[["totalmass"]]) && !is.null(Par[["mburst"]]) && !is.null(Par[["mSFR"]]) && is.null(Par[["burstfraction"]])){
      # Option C: totalmass, mSFR and mburst known; burstfraction=NULL
      cat("Case C\n") # ToDo: Improve verbosity
      # Nothing to touch on
    } else if (is.null(Par[["totalmass"]]) && !is.null(Par[["mburst"]]) && !is.null(Par[["mSFR"]]) && is.null(Par[["burstfraction"]])){
      # Option D: mSFR and mburst known; totalmass=burstfraction=NULL
      cat("Case D\n") # ToDo: Improve verbosity
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
  
  # set.seed(Parameters$RndSeed)
  
  a <- uniqueArguments(Parameters, speclib)
  
  
  
  
  
  
  
  
  
  return(a)
  
  
}

a <- generateDataFrameArguments(Parameters, n.simul=10, speclib=EMILESCombined)

# print(a)
agevec = EMILESCombined$Age
mass_val = splat(a$massfunc)(c(list(agevec), a[names(a) %in% names(formals(a$massfunc))]))
plot(agevec, mass_val, type="l", main=a$burstfraction)
abline(v=13.8e9, col="red")

print(a$totalmass)
print(a$mSFR)
print(a$mburst)
print(a$mburstage)
print(a$burstfraction)
cat("---\n")
print(a$mskew)
print(a$mperiod)




# 
# 
# agevec = EMILESCombined$Age
# ageweights=EMILESCombined$AgeWeights
# mf1 <- massfunc_snorm_burst(agevec, mSFR=1)*EMILESCombined$AgeWeights
# mf2 <- massfunc_snorm_burst(agevec, mSFR=1, mburst = 1e-4)*EMILESCombined$AgeWeights
# 
# sum(mf1)
# sum(mf2)
# sum(mf2) - sum(mf1)
# (sum(mf2)-sum(mf1))/sum(mf2)
# 
# 
# mburstage=0.1
# mburstage_scaled = mburstage*1e9
# mSFR1 = 1
# 
# # 1
# M_total_SFR1 = sum(massfunc_snorm_burst(agevec, mSFR1)*ageweights)
# print(M_total_SFR1)
# print(log10(M_total_SFR1))
# 
# # 2
# totalmass = 16521794      # Sacado de la distribución
# mburst_dist = 0.000445
# 
# 
# # 3
# M_burst = sum(mburst_dist * ageweights[agevec<mburstage_scaled])
# 
# fraction = M_burst/totalmass
# print(totalmass)
# print(log10(totalmass))
# mSFR_adjusted = totalmass*(1-fraction)/M_total_SFR1
# 
# 
# dev.off()
# mass_val = massfunc_snorm_burst(agevec, mSFR=mSFR_adjusted, mburst=mburst_dist, mburstage=mburstage)
# plot(agevec, mass_val, type="l", log="x")
# 
# plot(mass_val)


