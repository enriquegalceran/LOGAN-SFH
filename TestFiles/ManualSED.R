


setwd("~/Documents/GitHub/LOGAN-SFH")
source("FuncionesAuxiliaresYLibrerias.R")

filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filtersData <- list()
for (filter in filtersHST) {
  filtersData[[filter]] <-
    read.table(
      paste("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat", sep = ""),
      col.names = c("wave", "response")
    )
}

data("EMILES")
data("BC03hr")

##############################
if (TRUE){SFH = SFHfunc
  z = 0.1
  tau_birth = 1
  tau_screen = 0.3
  tau_AGN = 1
  pow_birth = -0.7
  pow_screen = -0.7
  pow_AGN = -0.7
  alpha_SF_birth = 1
  alpha_SF_screen = 3
  alpha_SF_AGN = 0
  AGNlum = 0
  sparse = 5
  Dale = NULL
  AGN = NULL
  filters = "all"
  Dale_M2L_func = NULL
  returnall = TRUE
  H0 = 67.8
  OmegaM = 0.308
  OmegaL = 1 - OmegaM
  waveout = seq(2, 9.35, by = 0.01)
  unimax = 1.38e+10
  agemax = NULL
  LumDist_Mpc = NULL
  addradio = FALSE
  Te = 10000
  ff_frac = 0.1
  ff_power = -0.1
  sy_power = -0.8
  AGNct = 60
  AGNal = 4
  AGNbe = -0.5
  AGNta = 1
  AGNrm = 60
  AGNan = 30
  Eb = 0
  L0 = 2175.8
  LFWHM = 470
  IGMabsorb = 0
}##############################


{
  massfunc = massfunc_dtau
  mSFR = 100
  mpeak = 10
  tau_birth = 1.0
  z = 0.0
  tau_screen = 0.3
  filtout = filtersData
  sparse = 1
  speclib = EMILES
  emission = TRUE
  emission_scale = 'SFR'
  
  EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
  speclib1=EMILES
  speclib2=EMILESCombined
}

# SFH y espectros de momento
Stars = SFHfunc(
  z = z,
  tau_birth = tau_birth,
  tau_screen = tau_screen,
  pow_birth = pow_birth,
  pow_screen = pow_screen,
  sparse = sparse,
  speclib = speclib1,
  filters = filtersData,
  unimax = unimax,
  agemax = agemax,
  Eb = Eb,
  L0 = L0,
  LFWHM = LFWHM,
  massfunc = massfunc_dtau,
  mSFR = 100,
  mpeak = 10,
  emission = TRUE,
  emission_scale = emission_scale,
  Z=0.02
)

Stars2 = SFHfunc(
  z = z,
  tau_birth = tau_birth,
  tau_screen = tau_screen,
  pow_birth = pow_birth,
  pow_screen = pow_screen,
  sparse = sparse,
  speclib = speclib2,
  filters = filtersData,
  unimax = unimax,
  agemax = agemax,
  Eb = Eb,
  L0 = L0,
  LFWHM = LFWHM,
  massfunc = massfunc_dtau,
  mSFR = 100,
  mpeak = 10,
  emission = TRUE,
  emission_scale = emission_scale,
  Z=0.02
)


plot(Stars$wave_lum, Stars$lum_unatten, log='xy', xlim=c(1.5e3,6e4), ylim=c(1e5,1e8),
     xlab='Wavelength (Ang)', ylab='Lum (Lsol/Ang)', type="l", col="blue",
     main="Atten/Unatten-EMILES(Old/Young+Old) - emissions ON")
lines(Stars$wave_lum, Stars$lum_atten, col='green')
lines(Stars2$wave_lum, Stars2$lum_unatten, col='red')
lines(Stars2$wave_lum, Stars2$lum_atten, col='maroon')
abline(v=wavelengthRangeInMUSE[1], col="black")
abline(v=wavelengthRangeInMUSE[2], col="black")
legend('topright',
       legend=c("Unatt-O", "Atten-O", "Unatt-Y+O", "Atten-Y+O", "MUSE"),
       col=c('blue', 'green', 'red', 'maroon', 'black'),
       lty=c(1, 1, 1, 1), cex=0.75
)



# SFH y espectros de momento
Stars3 = SFHfunc(
  z = z,
  tau_birth = tau_birth,
  tau_screen = tau_screen,
  pow_birth = pow_birth,
  pow_screen = pow_screen,
  sparse = sparse,
  speclib = speclib1,
  filters = filtersData,
  unimax = unimax,
  agemax = agemax,
  Eb = Eb,
  L0 = L0,
  LFWHM = LFWHM,
  massfunc = massfunc_dtau,
  mSFR = 100,
  mpeak = 10,
  emission = FALSE,
  emission_scale = emission_scale,
  Z=0.02
)

Stars4 = SFHfunc(
  z = z,
  tau_birth = tau_birth,
  tau_screen = tau_screen,
  pow_birth = pow_birth,
  pow_screen = pow_screen,
  sparse = sparse,
  speclib = speclib2,
  filters = filtersData,
  unimax = unimax,
  agemax = agemax,
  Eb = Eb,
  L0 = L0,
  LFWHM = LFWHM,
  massfunc = massfunc_dtau,
  mSFR = 100,
  mpeak = 10,
  emission = FALSE,
  emission_scale = emission_scale,
  Z=0.02
)


plot(Stars3$wave_lum, Stars3$lum_unatten, log='xy', xlim=c(1.5e3,6e4), ylim=c(1e5,1e8),
     xlab='Wavelength (Ang)', ylab='Lum (Lsol/Ang)', type="l", col="blue",
     main="Atten/Unatten-EMILES(Old/Young+Old) - emissions OFF")
lines(Stars3$wave_lum, Stars3$lum_atten, col='green')
lines(Stars4$wave_lum, Stars4$lum_unatten, col='red')
lines(Stars4$wave_lum, Stars4$lum_atten, col='maroon')
abline(v=wavelengthRangeInMUSE[1], col="black")
abline(v=wavelengthRangeInMUSE[2], col="black")
legend('topright',
       legend=c("Unatt-O", "Atten-O", "Unatt-Y+O", "Atten-Y+O", "MUSE"),
       col=c('blue', 'green', 'red', 'maroon', 'black'),
       lty=c(1, 1, 1, 1), cex=0.75
)


plot(Stars$agevec, Stars$Zvec, main="metal Stars")
plot(Stars2$agevec, Stars2$Zvec, main="metal Stars2")
plot(Stars3$agevec, Stars3$Zvec, main="metal Stars3")
plot(Stars4$agevec, Stars4$Zvec, main="metal Stars4")



##########################


if (TRUE){
  stellpop = EMILES
  namestellpop = "EMILES"
  colvec = rev(rainbow(length(seq(1, length(stellpop$Age), by=2)), end = 2 / 3))
  
  # for (metal in 1:length(stellpop$Z)){
  for (metal in 1:1){
    fprintf("Metal: %f\n", stellpop$Z[metal])
    ymin = min(stellpop$Zspec[[metal]][stellpop$Zspec[[metal]] != 0])
    ymax = max(stellpop$Zspec[[metal]][stellpop$Zspec[[metal]] != 0])
    
    plot(stellpop$Wave, stellpop$Zspec[[metal]][1,], log="xy", type="l", col=colvec[1], ylim=c(ymin, ymax),
         main=sprintf("Metallicity: %f. (Blue=Youngest, Red=Oldest)", EMILES$Z[metal]), xlab='Wavelength (Ang)', ylab='Response')
    for (age in seq(3, length(stellpop$Age), by=2)){
      lines(stellpop$Wave, stellpop$Zspec[[metal]][age,], col=colvec[age])
    }
    abline(v=wavelengthRangeInMUSE[1], col="black")
    abline(v=wavelengthRangeInMUSE[2], col="black")
  }
}


if (TRUE){
  xLim = c(min(c(EMILES$Z, EMILESCombined$Z)), max(c(EMILES$Z, EMILESCombined$Z)))
  yLim = c(min(c(EMILES$Age, EMILESCombined$Age)), max(c(EMILES$Age, EMILESCombined$Age)))
  
  plot(rep(EMILES$Z[1], length(EMILES$Age)), EMILES$Age,
       log="y", ylab = "Age / Yrs", xlab="Metallicity", xlim=xLim, ylim=yLim, pch=16, col="red",
       main='"Resolution" of data')
  for (Zi in 1:length(EMILES$Z)){
    points(rep(EMILES$Z[Zi], length(EMILES$Age)), EMILES$Age, pch=16, col="red")
  }
  
  points(rep(EMILESCombined$Z[1], length(EMILESCombined$Age)), EMILESCombined$Age,
         ylab = "Age / Yrs", xlab="Metallicity", pch=16, col="blue")
  for (Zi in 1:length(EMILESCombined$Z)){
    points(rep(EMILESCombined$Z[Zi], length(EMILESCombined$Age)), EMILESCombined$Age, pch=16, col="blue")
  }
  
  legend("bottomright", legend=c("EMILES", "EMILESCombined"), pch=16, col=c("red", "blue"))
}

if (TRUE){
  CommonZ = c(0.0001, 0.001, 0.004, 0.030)
  CommonAgeSum = c()
  wiggleroom = 0.02
  for (age in EMILESCombined$Age){
    CommonAgeSum = c(CommonAgeSum, sum(EMILES$Age > age*(1-wiggleroom) & EMILES$Age < age*(1+wiggleroom)))
  }
  CommonAgeE = c()
  CommonAgeEC = c()
  for (age in 1:length(EMILESCombined$Age)){
    if (CommonAgeSum[age] == 1){
      CommonAgeE = c(CommonAgeE, EMILES$Age[EMILES$Age > EMILESCombined$Age[age]*(1-wiggleroom) & EMILES$Age < EMILESCombined$Age[age]*(1+wiggleroom)])
      CommonAgeEC = c(CommonAgeEC, EMILESCombined$Age[age])
    }
  }
  CommonAgeEC = CommonAgeEC[-33]      # No entiendo por que aparece uno adicional, no deberia...
  CommonAgeE - CommonAgeEC
  CommonAgeE = CommonAgeE[6, 16, 29]
  CommonAgeEC = CommonAgeEC[6, 16, 29]
  # CommonZE = c(0.0001, 0.003, 0.001, 0.004, 0.0198, 0.030)
  # CommonZEC = c(0.0001, 0.004, 0.001, 0.004, 0.0190, 0.030)
  CommonZE = c(0.0003)
  CommonZEC = c(0.0004)
  
  for (agei in 1:length(CommonAgeE)){
    ageE = CommonAgeE[agei]
    ageEC = CommonAgeEC[agei] 
    
    colvecE = rev(rainbow(length(seq(1, length(CommonZE))), end = 1 / 3))
    colvecEC = rev(rainbow(length(seq(1, length(CommonZEC))), start=1/2, end = 1 / 3))
    
    
    ylimData = c()
    for (zi in 1:length(CommonZE)){
      ylimDatatmp = c()
      ylimDatatmp = c(ylimDatatmp, EMILES$Zspec[[match(CommonZE[zi], EMILES$Z)]][match(ageE, EMILES$Age),])
      ylimDatatmp = c(ylimDatatmp, EMILESCombined$Zspec[[match(CommonZEC[zi], EMILESCombined$Z)]][match(ageEC, EMILESCombined$Age),])
      ylimData = c(ylimData, ylimDatatmp[ylimDatatmp != 0 & !is.na(ylimDatatmp)])
    }
    yLim = c(min(ylimData), max(ylimData))
    
    legendText = c()
    legendColor = c()
    for (zi in 1:length(CommonZE)){
      zE = CommonZE[zi]
      zEC = CommonZEC[zi]
      if (zi == 1){
        plot(EMILES$Wave, EMILES$Zspec[[match(zE, EMILES$Z)]][match(ageE, EMILES$Age),],
             log="xy", type="l", ylim=yLim, col=colvecE[zi],
             main=sprintf("Age: %4.2e", ageE))
        lines(EMILESCombined$Wave, EMILESCombined$Zspec[[match(zEC, EMILESCombined$Z)]][match(ageE, EMILESCombined$Age),], col=colvecEC[zi])
      } else {
        lines(EMILES$Wave, EMILES$Zspec[[match(zE, EMILES$Z)]][match(ageE, EMILES$Age),], col=colvecE[zi])
        lines(EMILESCombined$Wave, EMILESCombined$Zspec[[match(zEC, EMILESCombined$Z)]][match(ageE, EMILESCombined$Age),], col=colvecEC[zi])

      }
      legendText = c(legendText, sprintf("%6.4f-E", zE), sprintf("%6.4f-EC", zEC))
      legendColor = c(legendColor, colvecE[zi], colvecEC[zi])
    }
    
    legend("bottomright", legend=legendText, col=legendColor, lty=1)
    
  }
  
  
}


Zsol = 0.0198

MH= c(-2.32, -1.71, -1.31, -0.71, -0.4, 0, 0.22)

zmh = c()

for (mh in MH){
  zmh <- c(zmh, Zsol*10^mh)
}
zmh

EMILES2 = EMILES
EMILES2[["Z"]] = 0.0004



####################################################################################
# plot(Stars$wave_lum, Stars$lum_unatten, log='xy', xlim=c(1.5e3,6e4), ylim=c(1e5,1e8),
#      xlab='Wavelength (Ang)', ylab='Lum (Lsol/Ang)', type="l", col="blue",
#      main="Atten/Unatten-EMILES(Old) - emissions ON")
# lines(Stars$wave_lum, Stars$lum_atten, col='green')
# lines(Stars3$wave_lum, Stars3$lum_unatten, col='red')
# lines(Stars3$wave_lum, Stars3$lum_atten, col='maroon')
# abline(v=wavelengthRangeInMUSE[1], col="black")
# abline(v=wavelengthRangeInMUSE[2], col="black")
# legend('topright',
#        legend=c("Unatt-ON", "Atten-ON", "Unatt-OFF", "Atten-OFF", "MUSE"),
#        col=c('blue', 'green', 'red', 'maroon', 'black'),
#        lty=c(1, 1, 1, 1), cex=0.75
# )

# lines(testSED3$DustEmit, col='brown')
# lines(testSED3$AGN, col='purple')
# legend('topright',
#        legend=c('Total Lum', 'Star Un-Atten', 'Stars Atten', 'Dust Emit', 'AGN'),
#        col=c('black', 'blue', 'darkgreen', 'brown', 'purple'),
#        lty=c(1,2,1,1,1),
#        lwd=c(5,1,1,1,1)
# )




testObject = customWaveout(Stars, seq(5000, 10000))

limitsSpecLib = c(min(speclib$Wave), max(speclib$Wave))
insideRangeIdx = Stars$wave_lum >= limitsSpecLib[1] & Stars$wave_lum <= limitsSpecLib[2]
sum(insideRangeIdx)

sum(is.na(Stars$lum_unatten))
sum(is.na(Stars$lum_unatten[insideRangeIdx]))
lumtot_atten2 = sum(qdiffCustom(Stars$wave_lum[insideRangeIdx])*Stars$lum_atten[insideRangeIdx])
lumtot_atten2
Stars$lumtot_unatten

noNALumIdx = !is.na(Stars$lum_unatten)
sum(noNALumIdx)
sum(noNALumIdx) - length(noNALumIdx)

plot(Stars$wave_lum, Stars$lum_unatten, log='xy', xlim=c(1.5e3,6e4), ylim=c(1e5,1e8),
     xlab='Wavelength (Ang)', ylab='Lum (Lsol/Ang)', type="l", col="blue",
     main="Atten/Unatten-EMILES(Old/Young+Old) - emissions ON")
lines(Stars$wave_lum, Stars$lum_atten, col='green')
abline(v=wavelengthRangeInMUSE[1], col="black")
abline(v=wavelengthRangeInMUSE[2], col="black")

points(Stars$wave_lum[!noNALumIdx], rep(5e7, length(Stars$wave_lum[!noNALumIdx])))

legend('topright',
       legend=c("Unatt-O", "Atten-O", "Unatt-Y+O", "Atten-Y+O", "MUSE"),
       col=c('blue', 'green', 'red', 'maroon', 'black'),
       lty=c(1, 1, 1, 1), cex=0.75
)

plot(Stars$wave_lum, Stars$lum_unatten-Stars$lum_atten, log="xy", type="l")


Stars$wave_lum = Stars$wave_lum[insideRangeIdx]
Stars$lum_unatten = Stars$lum_unatten[insideRangeIdx]
Stars$lum_atten = Stars$lum_atten[insideRangeIdx]
Stars$lumtot_atten = sum(qdiffCustom(Stars$wave_lum) * Stars$lum_atten)



Dust_Birth = Dale_interp(alpha_SF = alpha_SF_birth, Dale = Dale)
Dust_Screen = Dale_interp(alpha_SF = alpha_SF_screen, Dale = Dale)
SED_Bdust_Sdust = Dust_Birth$Aspec * Stars$lumtot_birth +
  Dust_Screen$Aspec * Stars$lumtot_screen
waveDust = Dust_Birth$Wave
plot(waveDust, abs(SED_Bdust_Sdust), col='brown', log="xy", type="l")



SED_Stars_Bdust_Sdust = addspec(
  wave1 = Stars$wave_lum,
  flux1 = Stars$lum_atten,
  wave2 = Dust_Screen$Wave,
  flux2 = SED_Bdust_Sdust,
  extrap = 0,
  waveout = waveout
)


if (!is.null(Dale_M2L_func) & returnall) {
  dustlum_birth = Stars$lumtot_birth
  dustlum_screen = Stars$lumtot_screen
  dustmass_birth = Stars$lumtot_birth / Dale_M2L_func(alpha_SF_birth)
  dustmass_screen = Stars$lumtot_screen / Dale_M2L_func(alpha_SF_screen)
} else {
  dustlum_birth = 0
  dustlum_screen = 0
  dustmass_birth = 0
  dustmass_screen = 0
}
if (is.null(AGN) | AGNlum == 0) {
  Final = SED_Stars_Bdust_Sdust
  AGN = NULL
  dustlum_AGN = 0
  dustmass_AGN = 0
} else {
  if (inherits(AGN, "Fritz")) {
    AGN = AGNinterp(
      lum = AGNlum,
      ct = AGNct,
      al = AGNal,
      be = AGNbe,
      ta = AGNta,
      rm = AGNrm,
      an = AGNan,
      Fritz = AGN
    )
    dustlum_AGN = NA
    dustmass_AGN = NA
    AGN = atten_emit(
      wave = AGN$wave,
      flux = AGN$lum *
        .erg_to_lsol,
      tau = tau_screen,
      pow = pow_screen,
      alpha_SF = alpha_SF_screen,
      Dale = Dale,
      Dale_M2L_func = Dale_M2L_func,
      waveout = waveout,
      Eb = Eb,
      L0 = L0,
      LFWHM = LFWHM
    )
    if (!is.null(Dale_M2L_func) & returnall) {
      dustlum_screen = dustlum_screen + AGN$total_atten
      dustmass_screen = dustmass_screen + AGN$dustmass
    }
    AGN = AGN$final
    Final = data.frame(wave = SED_Stars_Bdust_Sdust$wave,
                       flux = SED_Stars_Bdust_Sdust$flux + AGN$flux)
    colnames(AGN)[2] = "lum"
  } else {
    AGN = atten_emit(
      wave = AGN$Wave,
      flux = AGN$Aspec *
        AGNlum * .erg_to_lsol,
      tau = tau_AGN,
      pow = pow_AGN,
      alpha_SF = alpha_SF_AGN,
      Dale = Dale,
      Dale_M2L_func = Dale_M2L_func,
      waveout = waveout
    )
    if (!is.null(Dale_M2L_func) & returnall) {
      dustlum_AGN = AGN$total_atten
      dustmass_AGN = AGN$dustmass
    }
    AGN = atten_emit(
      wave = AGN$final$wave,
      flux = AGN$final$flux,
      tau = tau_screen,
      pow = pow_screen,
      alpha_SF = alpha_SF_screen,
      Dale = Dale,
      Dale_M2L_func = Dale_M2L_func,
      waveout = waveout,
      Eb = Eb,
      L0 = L0,
      LFWHM = LFWHM
    )
    if (!is.null(Dale_M2L_func) & returnall) {
      dustlum_screen = dustlum_screen + AGN$total_atten
      dustmass_screen = dustmass_screen + AGN$dustmass
    } else {
      dustlum_AGN = 0
      dustmass_AGN = 0
    }
    AGN = AGN$final
    if (length(SED_Stars_Bdust_Sdust$flux) == length(AGN$flux)) {
      Final = data.frame(wave = SED_Stars_Bdust_Sdust$wave,
                         flux = SED_Stars_Bdust_Sdust$flux + AGN$flux)
    } else {
      Final = addspec(
        wave1 = SED_Stars_Bdust_Sdust$wave,
        flux1 = SED_Stars_Bdust_Sdust$flux,
        wave2 = AGN$wave,
        flux2 = AGN$flux
      )
    }
    colnames(AGN)[2] = "lum"
  }
}
if (addradio) {
  Final = radiocont(
    wave = Final$wave,
    flux = Final$flux,
    z = 0,
    Te = Te,
    ff_frac = ff_frac,
    ff_power = ff_power,
    sy_power = sy_power,
    flux_in = "wave",
    flux_out = "wave"
  )
}
colnames(Final)[2] = "lum"
if (IGMabsorb > 0) {
  sel = which(Final$wave < 1215.67)
  Final$lum[sel] = Final$lum[sel] * (1 - IGMabsorb)
  sel = which(Final$wave < 911.8)
  Final$lum[sel] = 0
}


lines(Final, col="black")


if (is.null(filtout) & !is.null(filters)) {
  if (filters[1] == "all") {
    cenwave = NULL
    data("cenwave", envir = environment())
    filters = cenwave$filter
  }
  filtout = list()
  for (i in filters) {
    filtout = c(filtout, list(getfilt(i)))
  }
}
if (z > 0 & !is.null(filtout)) {
  Flux = Lum2Flux(
    wave = Final$wave,
    lum = Final$lum,
    z = z,
    H0 = H0,
    OmegaM = OmegaM,
    OmegaL = OmegaL,
    ref = ref,
    LumDist_Mpc = LumDist_Mpc
  )
  Flux$flux = convert_wave2freq(flux_wave = Flux$flux *
                                  .cgs_to_jansky,
                                wave = Flux$wave)
  photom_out = {
    
  }
  for (i in 1:length(filtout)) {
    photom_out = c(photom_out,
                   bandpass(
                     flux = Flux$flux,
                     wave = Flux$wave,
                     filter = filtout[[i]]
                   ))
  }
} else if (z > 0 & is.null(filtout)) {
  Flux = Lum2Flux(
    wave = Final$wave,
    lum = Final$lum,
    z = z,
    H0 = H0,
    OmegaM = OmegaM,
    OmegaL = OmegaL,
    ref = ref,
    LumDist_Mpc = LumDist_Mpc
  )
  Flux$flux = convert_wave2freq(flux_wave = Flux$flux *
                                  .cgs_to_jansky,
                                wave = Flux$wave)
  photom_out = NULL
} else if (z <= 0 & !is.null(filtout)) {
  Flux = cbind(wave = Final$wave,
               flux = Final$lum * .lsol_to_absolute)
  photom_out = photom_flux(Flux, outtype = "magAB", filters = filtout)
} else {
  Flux = NULL
  photom_out = NULL
}
if (returnall) {
  StarsAtten = data.frame(wave = Stars$wave_lum, lum = Stars$lum_atten)
  StarsUnAtten = data.frame(wave = Stars$wave, lum = Stars$lum_unatten)
  DustEmit = data.frame(wave = Dust_Screen$Wave, lum = SED_Bdust_Sdust)
  output = list(
    Photom = photom_out,
    FinalFlux = Flux,
    FinalLum = Final,
    StarsAtten = StarsAtten,
    StarsUnAtten = StarsUnAtten,
    DustEmit = DustEmit,
    AGN = AGN,
    Stars = Stars,
    dustmass = c(
      birth = dustmass_birth,
      screen = dustmass_screen,
      AGN = dustmass_AGN,
      total = sum(c(
        dustmass_birth, dustmass_screen,
        dustmass_AGN
      ), na.rm = TRUE)
    ),
    dustlum = c(
      birth = dustlum_birth,
      screen = dustlum_screen,
      AGN = dustlum_AGN,
      total = sum(c(
        dustlum_birth, dustlum_screen,
        dustlum_AGN
      ), na.rm = TRUE)
    ),
    call = call
  )
  class(output) = "ProSpectSED"
  return(output)
} else {
  return(photom_out)
}
