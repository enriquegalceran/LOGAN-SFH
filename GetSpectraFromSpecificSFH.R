#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

original_parameters <- par()
library(ProSpect)
library(plotrix)
source("MUTANTS/LOGAN-CLAWS.R")
setwd("~/Documents/GitHub/LOGAN-SFH")
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")

# csvData <- read.csv(file = args[[0]])
csvData <- read.csv(file = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file.pd")
csvData2 <- read.csv(file = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file2.pd")
csvData3 <- read.csv(file = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file3.pd")
csvData4 <- read.csv(file = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file4.pd")
waveout=seq(4700, 9400, 1.25)
filt_cenwave = c(2708.680, 3358.580, 4329.475, 5331.898, 8072.917)
# col_magnitudes = col2rgb(w_length2rgb(filt_cenwave/10))

configFilename = "Data_Generation_Parameters.R"
source(configFilename)


#### Plot input data ####
plot(csvData$sfh_true)
plot(csvData$z_true)
plot(csvData$agevec, csvData$z0, log="x", type="l")


simple.two.axes.plot <- function(lx, ly, rx=lx, ry, xlab=NULL, lylab=NULL, rylab=NULL, main="",
                                 ltype="l", rtype="p", lpch=1, rpch=17, lcol="black", rcol="red", xlim=NULL,
                                 log=""){
    if (is.null(xlim)){
        xlim = c(min(c(lx, rx)), max(c(lx, rx)))
    }
    par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
    plot(lx, ly, type=ltype, xlim=xlim,
         main=main, xlab=xlab, ylab=lylab, pch=lpch, col=lcol, log=log) # first plot
    par(new = TRUE)
    plot(rx, ry, type=rtype, pch=rpch, axes = FALSE, bty = "n", xlab = "", ylab = "", col=rcol, xlim=xlim, log=log)
    axis(side=4, at = pretty(range(ry)), col=rcol, col.axis=rcol)
    mtext("magnitudes", side=4, line=3, col=rcol)
}


simple.two.axes.plot(waveout, csvData2$spectr_in, filt_cenwave, csvData3$magnitudes_in,
                     xlab="Wavelength (A)", lylab="Flux", rylab="Magnitudes",
                     main="Spectra + Magnitudes - Input data")


filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filters <- list()
for (filter in filtersHST) {
    filters[[filter]] <-
        read.table(
            paste0("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat"),
            col.names = c("wave", "response")
        )
}


#### Define custom functions ####
massfunc_custom <- function (age, magevec, msfh, magemax = 13.8){
    magemax = magemax * 1e+09
    age[age < 1e+05] = 1e+05
    temp = splinefun(log10(magevec), msfh, method = 'monoH.FC')(log10(age))
    temp[temp < 0] = 0
    temp[age > magemax] = 0
    return(temp)
}

Z_custom <- function (age, Zagevec, Zsfh, Zagemax = 13.8, ...){
    Zagemax = Zagemax * 1e+09
    age[age < 1e+05] = 1e+05
    temp = splinefun(log10(Zagevec), Zsfh, method = 'monoH.FC')(log10(age))
    temp[temp < 0] = 0
    temp[age > Zagemax] = 0
    return(temp)
}

post_process_spectra <- function(object, waveout){
    # SFR
    tmp_output = convertAgevecToOutputScale(object$agevec, object$SFR)
    object$SFR = tmp_output$data
    # massvec
    tmp_output = convertAgevecToOutputScale(object$agevec, object$massvec)
    object$massvec = tmp_output$data
    # Zvec
    tmp_output = convertAgevecToOutputScale(object$agevec, object$Zvec)
    object$Zvec = tmp_output$data
    # agevec
    tmp_output = convertAgevecToOutputScale(object$agevec, object$agevec)
    object$agevec = suppressWarnings(tmp_output$data)    
    
    # Adjust Spectra to Wavelength (waveout)
    object$flux = interpolateToWaveout(
        lapply(object$flux["wave"], as.numeric)[[1]],
        lapply(object$flux["flux"], as.numeric)[[1]],
        waveout,
        returnList=TRUE)
    
    return(object)
}

compare_two_objects <- function(object1, object2, main="Compare 2 objects"){
    par(mfrow=c(2,2))
    simple.two.axes.plot(
        object1$flux$wave, object1$flux$flux,
        object2$flux$wave, object2$flux$flux,
        main="Spectra", ltype="l", rtype="l",
        lylab="flux", rylab="flux", xlab="Wavelength (A)"
    )
    simple.two.axes.plot(
        object1$out$cenwave, object1$out$out,
        object2$out$cenwave, object2$out$out,
        main="Magnitudes", ltype="p", rtype="p", lpch=1, rpch=2,
        lylab="magnitudes", rylab="magnitudes", xlab="Wavelength (A)"
    )
    simple.two.axes.plot(
        object1$agevec, object1$SFR,
        object2$agevec, object2$SFR,
        main="SFR", ltype="l", rtype="l", log="x",
        lylab="SFR", rylab="SFR", xlab="Lookbacktime (yrs)"
    )
    simple.two.axes.plot(
        object1$agevec, object1$Zvec,
        object2$agevec, object2$Zvec,
        main="Z", ltype="l", rtype="l", log="x",
        lylab="Z", rylab="Z", xlab="Lookbacktime (yrs)"
    )
    par(mfrow=c(1,1))
    mtext(main, side=3, line=-1, outer=TRUE)
}


#### Calculate spectra ####
{
    spectra_no_stand = SFHfunc(
        massfunc = massfunc_custom,
        speclib = EMILESCombined,
        filters=filters,
        Z = Z_custom,
        magevec=csvData[["agevec"]],
        msfh=csvData[["sfh_no_stand"]],
        Zagevec = csvData[["agevec"]],
        Zsfh=csvData[["z_no_stand"]],
    )
    spectra_no_stand <- post_process_spectra(spectra_no_stand, waveout)
    
    spectra_true = SFHfunc(
        massfunc = massfunc_custom,
        speclib = EMILESCombined,
        filters=filters,
        Z = Z_custom,
        magevec=csvData[["agevec"]],
        msfh=csvData[["sfh_true"]],
        Zagevec = csvData[["agevec"]],
        Zsfh=csvData[["z_true"]],
    )
    spectra_true <- post_process_spectra(spectra_true, waveout)
    
    compare_two_objects(spectra_true, spectra_no_stand,
                        main=paste0("true vs no_stand - ID=", csvData3$ID[1],
                                    " - ", csvData4[[paste0("name", idx)]]))
    
    for (idx in 0:3){
        spectraObject = SFHfunc(
            massfunc = massfunc_custom,
            speclib = EMILESCombined,
            filters=filters,
            Z = Z_custom,
            magevec=csvData[["agevec"]],
            msfh=csvData[[paste0("sfh", idx)]],
            Zagevec = csvData[["agevec"]],
            Zsfh=csvData[[paste0("z", idx)]],
        )
        spectraObject <- post_process_spectra(spectraObject, waveout)
        
        
        compare_two_objects(spectra_true,
                            spectraObject,
                            main=paste0("true vs idx=", idx,
                                        " - ID=", csvData3$ID[1],
                                        " - ", csvData4[[paste0("name", idx)]]))
    }
}



simple.two.axes.plot(spectraObject$agevec, spectraObject$Zvec,
                     csvData[["agevec"]], csvData[[paste0("z", idx)]],
                     main="Z spectraObject vs SFH data",
                     ltype="l", rtype="l", log="x")

simple.two.axes.plot(spectraObject$flux$wave, spectraObject$flux$flux,
                     spectraObject$out$cenwave, spectraObject$out$out,
                     main="spectraObject")
simple.two.axes.plot(waveout, csvData2$spectr_in,
                     filt_cenwave, csvData3$magnitudes_in,
                     main="Input to CNN")

simple.two.axes.plot(spectraObject$flux$wave, spectraObject$flux$flux)







