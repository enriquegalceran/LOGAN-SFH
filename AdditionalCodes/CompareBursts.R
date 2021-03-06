# Compare spectra for bursts

# Initialize variables
.pardefault <- par()
setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
library(RColorBrewer)
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
'%!in%' <- function(x,y)!('%in%'(x,y))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filters <- list()
for (filter in filtersHST) {
  filters[[filter]] <-
    read.table(
      paste("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat", sep = ""),
      col.names = c("wave", "response")
    )
}


pathoutputs="Plots/ComparisonBetweenCompleteSFHAndOnlyBursts/"
massfuncprefix="f_"
parameters_function = list(massfunc=massfunc_snorm_burst, mperiod=1, mpeak=12, mSFR=10, mburstage=0.1)
mburst_value=3
age = EMILESCombined$Age
wavelengthRangeInMUSE <- c(4750, 9350) 
x.lim = c(6.3e+06, 13.8e+09)
y.lim = c(0,parameters_function[["mSFR"]])
log=""


#####################
## Calculate spectra
spectraObjectComplete <- do.call('SFHfunc', c(parameters_function,
                                              list(speclib=EMILESCombined,
                                                   emission=TRUE,
                                                   Z=0.02,
                                                   emission_scale="SFR",
                                                   filters=filters,
                                                   mburst=mburst_value
                                                   )
                                              )
                                 )

# Get the mass of the burst
massVectorsInBurst <- spectraObjectComplete$massvec[spectraObjectComplete$agevec<(parameters_function$mburstage*1e+09)]
massOfBurst <- sum(massVectorsInBurst)

spectraObjectReduced  <- do.call('SFHfunc', c(parameters_function,
                                              list(speclib=EMILESCombined,
                                                   emission=TRUE,
                                                   emission_scale="SFR",
                                                   Z=0.02,
                                                   forcemass=massOfBurst,
                                                   mburst=mburst_value,
                                                   magemax=parameters_function[["mburstage"]],
                                                   filters=filters
                                                   )
                                              )
                                 )
spectraObjectNoBurst  <- do.call('SFHfunc', c(parameters_function,
                                              list(speclib=EMILESCombined,
                                                   emission=TRUE,
                                                   emission_scale="SFR",
                                                   Z=0.02,
                                                   mburst=0,
                                                   filters=filters
                                                   )
                                              )
                                 )

#####################
## Plot SFR
pdf(file=paste0(pathoutputs, massfuncprefix, "SFHComplete", ".pdf"), width=15, height=9)
plot(spectraObjectComplete$agevec, spectraObjectComplete$SFR,
     type="l", log=log, main="SFH Complete", ylim=y.lim, xlim=x.lim, xlab="Age", ylab="SFR")
dev.off()
# y_theo = do.call(parameters_function[["massfunc"]],
#                  c(list(age=age), parameters_function[names(parameters_function) %in% names(formals(parameters_function[["massfunc"]]))]) )
# plot(age, y_theo, log=log, type="l", main="SFH Theo", ylim=y.lim, xlim=x.lim, xlab="Age", ylab="SFR")
pdf(file=paste0(pathoutputs, massfuncprefix, "SFHReduced", ".pdf"), width=15, height=9)
plot(spectraObjectReduced$agevec, spectraObjectReduced$SFR,
     type="l", log=log, main="SFH Reduced", ylim=y.lim, xlim=x.lim, xlab="Age", ylab="SFR")
dev.off()
pdf(file=paste0(pathoutputs, massfuncprefix, "SFHNoBurst", ".pdf"), width=15, height=9)
plot(spectraObjectNoBurst$agevec, spectraObjectNoBurst$SFR,
     type="l", log=log, main="SFH Reduced", ylim=y.lim, xlim=x.lim, xlab="Age", ylab="SFR")
dev.off()

#####################
## Plot SFR
plotSFRlog=""
muse.logi <- (spectraObjectComplete$wave_lum > wavelengthRangeInMUSE[1]) & (spectraObjectComplete$wave_lum < wavelengthRangeInMUSE[2])
maxy = max(c(spectraObjectComplete$lum_atten[muse.logi],
             spectraObjectComplete$lum_unatten[muse.logi],
             spectraObjectReduced$lum_atten[muse.logi],
             spectraObjectReduced$lum_unatten[muse.logi],
             spectraObjectNoBurst$lum_atten[muse.logi],
             spectraObjectNoBurst$lum_unatten[muse.logi]),
           na.rm=TRUE)
miny = min(c(spectraObjectComplete$lum_atten[muse.logi],
             spectraObjectComplete$lum_unatten[muse.logi],
             spectraObjectReduced$lum_atten[muse.logi],
             spectraObjectReduced$lum_unatten[muse.logi],
             spectraObjectNoBurst$lum_atten[muse.logi],
             spectraObjectNoBurst$lum_unatten[muse.logi]),
           na.rm=TRUE)
y.lim = c(miny, maxy)
pdf(file=paste0(pathoutputs, massfuncprefix, "LuminosityComparison-AttenvsUnatten", ".pdf"), width=15, height=9)
plot(spectraObjectComplete$wave_lum[muse.logi], spectraObjectComplete$lum_atten[muse.logi],
     type="l", col=col_vector[1], log=plotSFRlog, xlab="Wavelength (A)", ylab="Luminosity(Lsol / Angstrom)", ylim=y.lim,
     main="Luminosity Comparison - Attenuated vs Unattenuated")
lines(spectraObjectComplete$wave_lum[muse.logi], spectraObjectComplete$lum_unatten[muse.logi], type="l", col=col_vector[2], log=plotSFRlog)
legend(x="topright", legend=c("Lum Attenuated", "Lum Unattenuated"), col=c(col_vector[1], col_vector[2]), lty=1, lwd=1)
dev.off()

## Attenuated Comparison
pdf(file=paste0(pathoutputs, massfuncprefix, "LuminosityComparisonAtten-CompletevsReduced", ".pdf"), width=15, height=9)
plot(spectraObjectComplete$wave_lum[muse.logi], spectraObjectComplete$lum_atten[muse.logi],
     type="l", col=col_vector[1], log=plotSFRlog, xlab="Wavelength (A)", ylab="Luminosity(Lsol / Angstrom)", ylim=y.lim,
     main="Luminosity Comparison Attenuated - Complete vs Burst")
lines(spectraObjectComplete$wave_lum[muse.logi], spectraObjectReduced$lum_atten[muse.logi], type="l", col=col_vector[3], log=plotSFRlog)
lines(spectraObjectComplete$wave_lum[muse.logi], spectraObjectNoBurst$lum_atten[muse.logi], type="l", col=col_vector[6], log=plotSFRlog)
legend(x="topright", legend=c("Lum Attenuated Complete", "Lum Attenuated Burst", "Lum Attenuated NoBurst"), col=c(col_vector[1], col_vector[3], col_vector[6]), lty=1, lwd=1)
dev.off()

## Unattenuated Comparison
pdf(file=paste0(pathoutputs, massfuncprefix, "LuminosityComparisonUnatten-CompletevsReduced", ".pdf"), width=15, height=9)
plot(spectraObjectComplete$wave_lum[muse.logi], spectraObjectComplete$lum_unatten[muse.logi],
     type="l", col=col_vector[2], log=plotSFRlog, xlab="Wavelength (A)", ylab="Luminosity(Lsol / Angstrom)", ylim=y.lim,
     main="Luminosity Comparison Unattenuated - Complete vs Burst")
lines(spectraObjectComplete$wave_lum[muse.logi], spectraObjectReduced$lum_unatten[muse.logi], type="l", col=col_vector[5], log=plotSFRlog)
lines(spectraObjectComplete$wave_lum[muse.logi], spectraObjectNoBurst$lum_unatten[muse.logi], type="l", col=col_vector[7], log=plotSFRlog)
legend(x="topright", legend=c("Lum Unattenuated Complete", "Lum Unattenuated Burst", "Lum Unattenuated NoBurst"), col=c(col_vector[2], col_vector[5], col_vector[7]), lty=1, lwd=1)
dev.off()

## Colors Comparison
y.lim = c(min(c(spectraObjectComplete$out$out, spectraObjectReduced$out$out, spectraObjectNoBurst$out$out)), 
          max(c(spectraObjectComplete$out$out, spectraObjectReduced$out$out, spectraObjectNoBurst$out$out)))
pdf(file=paste0(pathoutputs, massfuncprefix, "ColorComparison-CompletevsReduced", ".pdf"), width=15, height=9)
plot(spectraObjectComplete$out$cenwave, spectraObjectComplete$out$out, col=col_vector[8:12], pch=15, main="Color Comparison", ylim=y.lim)
points(spectraObjectReduced$out$cenwave, spectraObjectReduced$out$out, col=col_vector[8:12], pch=17)
points(spectraObjectReduced$out$cenwave, spectraObjectNoBurst$out$out, col=col_vector[8:12], pch=18)
legend(x="topright", legend=c(filtersHST, "Complete", "Burst", "NoBurst"),
       pch=c(rep(16,5), 15, 17, 18), col=c(col_vector[8:12], "black", "black", "black"))
dev.off()
























# Whenever this function is used, with a new massfunc_original AND SFHfunc is going
# to be used, the formals need to be updated.
# massfunc_cut <- function(age, mcut, massfunc_original, ...){
#   dots = list(...)
#   massfunc_args = dots[names(dots) %in% names(formals(massfunc_original))]
#   
#   mass_out = do.call("massfunc_original", c(list(age), massfunc_args))
#   
#   for (ageIdx in 1:length(age)){
#     if (age[ageIdx] > mcut * 1e+09){
#       mass_out[ageIdx] = 0
#     }
#   }
#   return(mass_out)
# }
# func = massfunc_snorm_burst
# formals(massfunc_cut) <- c(formals(massfunc_cut), formals(func)[names(formals(func)) %!in% names(formals(massfunc_cut))])
# 
# 
# y <- func(age, mburst=5, mburstage = 0.1)
# plot(age, y, type="l", log="x")
# y2 <- massfunc_cut(age, mcut=0.1, massfunc_snorm_burst, mburst=4, mburstage = 0.1)
# plot(age, y2, type="l", log="x")
