# Main code/file
# Light Output Generated by AstroNomical objects

"
  REMINDER FOR THE FUTURE:
  
  IF DEBUGGING IS GIVING OUT ERRORS (i.e. going to compiler::tryCmpFun), the 
  issue is because it is trying to debug the compiler (for a faster execution
  in the future). To fix this, the function needs to be 'precompiled':
    - set the working directory
    - Execute the source of LOGAN
    - Load libraries
    - COMMENT EVERY LIBRARY LINE OUT
    - Now it should work fine
  
"
setwd("~/Documents/GitHub/LOGAN-SFH")
source("MUTANTS/LOGAN.R")
library(ProSpect)
library(ggplot2)


# massParams = list(
#   dtau=list(
#     name="dtau",
#     func=massfunc_dtau,
#     mSFR=10,
#     mpeak=seq(7,14.2,0.3),
#     mtau=seq(0.5,3.5, 0.2)
#   ),
#   snorm=list(
#     name="snorm",
#     func=massfunc_snorm,
#     mSFR=10,
#     mpeak=seq(7,14.2,0.3),
#     mperiod=seq(0.5,1.5, 0.2),
#     mskew=seq(-0.5, 1, 0.1)
#   )
# )

massParams2 = list(
  dtau=list(
    name="dtau",
    func=massfunc_dtau,
    mpeak=c(seq(1,12,0.2)),
    mtau=seq(1, 4, 0.2)
  )
)

massParams3 = list(
  dtau=list(
    name="dtau",
    func=massfunc_dtau,
    mpeak=c(8, 9),
    mtau=2
  )
)


# massParams4 = list(
#   dtau=list(
#     name="dtau",
#     func=massfunc_dtau,
#     mSFR=10,
#     mpeak=c(7,8,14),
#     mtau=c(2,3)
#   ),
#   snorm=list(
#     name="snorm",
#     func=massfunc_snorm,
#     mSFR=10,
#     mpeak=c(7,8,9,10),
#     mperiod=5,
#     mskew=c(0.5, 1, 1.5)
#   )
# )


ZParams = list(
  func=Zfunc_massmap_box,
  Zstart=1e-4,
  yield=0.03,
  Zfinal=seq(0.02, 0.08, 0.02)
)


# ZParams2 = list(
#   func=Zfunc_massmap_box,
#   Zstart=1e-4,
#   yield=0.03,
#   Zfinal=0.02
# )


EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")

outputFolder = "OutputsGenerated"
outputFolder = "/Volumes/Elements/Outputs"
outputFolder = "/Users/enrique/Documents/GitHub/LOGAN-SFH/OutputFolder"
absolutePath = TRUE

# Delete files in outputfolder
if (FALSE){
  # ls = list.files(file.path(getwd(), outputFolder))
  ls = list.files(outputFolder)
  for (file in ls){
    # file.remove(file.path(getwd(), outputFolder, file))
    file.remove(file.path(outputFolder, file))
  }
}

ls1 = sort(list.files(outputFolder))

## Probar interpolar
if (FALSE){
  EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
  waveout = seq(4700, 9400, 1.25)
  filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
  filters <- list()
  for (filter in filtersHST) {
    filters[[filter]] <-
      read.table(
        paste("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat", sep = ""),
        col.names = c("wave", "response")
      )
  }
  
  Stars = SFHfunc(
    speclib = EMILESCombined,
    filters = filters,
    massfunc = massfunc_dtau,
    mSFR = 100,
    mpeak = 10,
    emission = TRUE,
    emission_scale = "SFR",
    Z=Zfunc_massmap_box,
    Zstart=1e-4,
    Zfinal=0.1
  )
  plot(Stars$flux$wave, Stars$flux$flux, type="l", log="xy",
       xlim=c(4650, 9500), ylim=c(8e-16, 2.5e-15), main="original")
  newy = interpolateToWaveout(x1=Stars$flux$wave,
                              y1=Stars$flux$flux,
                              waveout=waveout,
                              offset=0.5,
                              n.points.integrate = 50,
                              method.to.evaluate = "mean")
  plot(waveout, newy, type="l", log="xy",
       xlim=c(4650, 9500), ylim=c(8e-16, 2.5e-15), main="nuevo")
  
}



# wave = seq(4700, 9400, 1.25)
# agevec = c(6300000,7900000,10000000,12600000,15800000,20000000,25100000,
#            31600000,39800000,50100000,63100000,70800000,79400000,89100000,
#            100000000,112200000,125900000,141300000,158500000,177800000,
#            199500000,223900000,251200000,281800000,316200000,354800000,
#            398100000,446700000,501200000,562300000,631000000,707900000,
#            794300000,891300000,1000000000,1122000000,1258900000,1412500000,
#            1584900000, 1778300000,1995300000,2238700000,2511900000,
#            2818400000,3162300000,3548100000,3981100000,4466800000,
#            5011900000,5623400000,6309600000,7079500000,7943300000,
#            8912500000,10000000000,11220200000,12589300000,14125400000,
#            15848900000,17782800000)
randomSamples = 5
out <- generateSpecFromParams(massParams = massParams3,
                       ZParams = ZParams,
                       folderPath = outputFolder,
                       absolutePath = absolutePath,
                       randomSamples = randomSamples,
                       speclib = EMILESCombined,
                       confirmation = TRUE,
                       verbose=1,
                       verboseSteps=randomSamples + 1,
                       filters="default",
                       cleanOutputFolder=FALSE,
                       bytesForPython=10e6,
                       singleOutput=TRUE)
t = out$time
save(t, file=file.path(outputFolder, paste0("time_", out$name, ".rda")))
# load(file.path(outputFolder, "time.rda"))

plot(t[, 1], t[, 2], type="l")
n_t <- dim(t)[1]

n_gt <- dim(t)[1]/(randomSamples + 1)
x <- seq(1, n_gt)
gt <- numeric(n_gt)
e_gt_u <- numeric(n_gt)
e_gt_d <- numeric(n_gt)
i <- 1
for (j in seq(1, n_t - randomSamples, randomSamples + 1)){
  gt[i] <- mean(t[j:(j+randomSamples), 2])
  e_gt_u[i] <- t[j, 2]
  e_gt_d[i] <- mean(t[(j+1):(j+randomSamples-1), 2])
  i <- i + 1
}

m_rolling <- 10
rolling_average_t <- numeric(n_gt - m_rolling)
x_r <- seq(1:(n_gt - m_rolling))
for (i in 1:(n_gt - m_rolling)){
  rolling_average_t[i] <- mean(gt[i:(i+m_rolling)])
}

plot (x, gt, ylim=c(min(e_gt_d)-0.1,max(e_gt_u) + 0.1), log="y")
segments(x, e_gt_d, x, e_gt_u)
epsilon <- 0.1
segments(x-epsilon, e_gt_d, x+epsilon, e_gt_d)
segments(x-epsilon, e_gt_u, x+epsilon, e_gt_u)
lines(x_r, rolling_average_t)

plot(x, gt, type="line")
plot(t[,1], t[, 2], type="line")

data <- data.frame(x_r, rolling_average_t)
sp <- qplot(x, gt)+geom_errorbar(aes(x=x, ymin=e_gt_d, ymax=e_gt_u), width=0.25)
sp
# ToDo: How do you add a second line in a ggplot????¿?¿?¿?¿


# 
# # Verify total file size
# ls2 = sort(list.files(outputFolder))
# ls2 = ls2[!ls2 %in% ls1]
# 
# sizeI = 0
# sizeL = 0
# 
# for (file in ls2){
#   if (substr(file, 1, 5) == "Input"){
#     sizeI = sizeI + file.info(file.path(outputFolder, file))$size
#   } else {
#     sizeL = sizeL + file.info(file.path(outputFolder, file))$size
#   }
# }
# 
# roundMultiple <- function(number, multiple, f=round){
#   return(f(number/multiple)*multiple)
# }
# 
# bytes2Human <- function(B, decimals=2){
#   betterUnits = roundMultiple(log10(B), 3, trunc)
#   B = B/(10^betterUnits)
#   units = c("", "K", "M", "G", "T", "P", "E", "Z", "Y")
#   units = units[betterUnits/3 + 1]
#   return(paste0(toString(round(B, decimals)), units, "B"))
# }
# 
# sizeIs = bytes2Human(sizeI)
# sizeLs = bytes2Human(sizeL)
# sizeTotal = sizeI + sizeL
# # TODO: This is no longer valid. Maybe redo it, maybe just remove it...
# cat("Total Size: ", bytes2Human(sizeTotal), "\n", sep="")
# cat("Total DISK Size: ", bytes2Human(length(ls2)*262144),
#     " (+", length(ls2)*262000/sizeTotal*100, "%) [OUTDATED! NEEDS TO BE REMOVED OR FIXED. LOW PRIORITY]\n", sep="")
# # ToDo fix this section where it calculates disk space used. Recalculate disk space at the beginning as well.
# 
# dfRead <- readFITS(file.path(outputFolder, ls2[1]))






