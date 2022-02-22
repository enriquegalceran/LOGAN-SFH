# Setup LOGAN

# Windows:
setwd("~/GitHub/LOGAN-SFH")

# CONFIG
remove.small.files = FALSE


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
cat("Saving new file ...\n")
saveRDS(EMILESReconstr, "EMILESData/EMILESCombined.rds")

# Remove small files
if (remove.small.files){
  cat("Small files for EMILES are being removed ...\n")
  for (i in c(1:6, 8)){
    file.remove(paste0("EMILESData/EMILES", i, ".rds"))
  }
  for (j in 1:7){
    file.remove(paste0("EMILESData/EMILES7_", j, ".rds"))
  }
} else {
  warning("Small files for EMILES were not removed, but will not be used. It is adviced to remove them.")
}

cat("It is advised to set these files to unchanged, or cloning/commiting will mess this up...\n")
cat("1) Go to .../LOGAN-SFH\n")
cat("2) git update-index --assume-unchanged $(git ls-files | tr '\n' ' ')\n")
cat("3) cd EMILESDATA/\n")
cat("4) git update-index --assume-unchanged $(git ls-files | tr '\n' ' ')\n")
cat("5) Delete small files.")

