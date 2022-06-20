library(future.apply)
plan(multisession)

setwd("~/Documents/GitHub/LOGAN-SFH")

n <- 4L
iteration.times <- 3L

cat(paste0("Initiating Parallel processing with ", n, " threads...\n"))
for (i in 1:iteration.times){
  cat(paste0("Iteration #", i, " with ", n, " threads...\n"))
  empty <- future_replicate(n, {
      print(.Random.seed)
      source(file = "DataGeneration.R")
  })
}
