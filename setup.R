# Setup LOGAN

# CONFIG
remove.small.files = TRUE


# Read split EMILESData files
EMILESReconstr = list()
for (i in 1:6){
  EMILESReconstr = c(EMILESReconstr, readRDS(file=paste0("EMILESData/EMILES", i, ".rds")))
}
tmp = list()
for (j in 1:7){
  tmp = c(tmp, readRDS(file=paste0("EMILESData/EMILES7_", j, ".rds")))
}
EMILESReconstr = c(EMILESReconstr, list(Zspec=tmp))
EMILESReconstr = c(EMILESReconstr, readRDS(file="EMILESData/EMILES8.rds"))


# Save as a single file
saveRDS(EMILESCombined, "EMILESData/EMILESCombined.rds")

# Remove small files
if (remove.small.files){
  for (i in c(1:6, 8)){
    file.remove(paste0("EMILESData/EMILES", i, ".rds"))
  }
  for (j in 1:7){
    file.remove(paste0("EMILESData/EMILES7_", j, ".rds"))
  }
}

