library(future.apply)
plan(multisession)

setwd("~/Documents/GitHub/LOGAN-SFH")

n <- 3L

empty <- future_replicate(n, {
    print(.Random.seed)
    source(file = "DataGeneration.R")
})




# empty <- future_lapply(n.simulations, function(ii){
#     
#     source('DataGeneration.R')
#     
# }, future.seed = TRUE)