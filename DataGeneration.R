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


outputFolder = "/Volumes/Elements/TrainingData/"
configFilename = "Data_Generation_Parameters.R"

# Load Parameters from the Config File.
source(configFilename)





