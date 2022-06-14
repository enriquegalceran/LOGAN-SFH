library(future.apply)
plan(multisession)

setwd("~/Documents/GitHub/LOGAN-SFH")

n <- 4L

empty <- future_replicate(n, {
    print(.Random.seed)
    source(file = "DataGeneration.R")
})

