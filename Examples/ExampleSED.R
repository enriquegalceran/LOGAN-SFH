data('BC03lr')
data('Dale_NormTot')
data('AGN_UnOb_Sparse')
data('Dale_M2L_func')
filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA',
          'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE',
          'P100_Herschel', 'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')
filtout={}
for(i in filters){filtout=c(filtout,list(approxfun(getfilt(i))))}
#for(i in filters){filtout=c(filtout,list(getfilt(i)))} #This works too, but slower!

#Try playing around with the below to see what impact different AGN, dust etc have:

testSED=ProSpectSED(AGNlum=1e43, tau_birth=1, tau_screen=0.3, tau_AGN=2,
                    alpha_SF_birth=1, alpha_SF_screen=3, alpha_SF_AGN=0, speclib=BC03lr, Dale=Dale_NormTot,
                    AGN=AGN_UnOb_Sparse, filtout=filtersData, Dale_M2L_func=Dale_M2L_func, returnall=TRUE)

plot(testSED$FinalLum, log='xy', xlim=c(1e2,1e7), ylim=c(1e2,1e7),
     xlab='Wavelength (Ang)', ylab='Lum (Lsol/Ang)', type='l', lwd=5)
lines(testSED$StarsUnAtten, col='blue', lty=2)
lines(testSED$StarsAtten, col='green')
lines(testSED$DustEmit, col='brown')
lines(testSED$AGN, col='purple')
legend('topright',
       legend=c('Total Lum', 'Star Un-Atten', 'Stars Atten', 'Dust Emit', 'AGN'),
       col=c('black', 'blue', 'darkgreen', 'brown', 'purple'),
       lty=c(1,2,1,1,1),
       lwd=c(5,1,1,1,1)
)

# or easier:

plot(testSED)

#An example of a more physical SFH and ZH:

testSED2 = ProSpectSED(massfunc = massfunc_snorm, Z=Zfunc_massmap_box, filtout=filtersData)
plot(testSED2)

##################################
# This is new, not from the example
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
wavelengthRangeInMUSE <- c(4750, 9350)  # Ang
wavelengthRangeInMUSELog10 <- log10(wavelengthRangeInMUSE)
waveoutManual = seq(wavelengthRangeInMUSELog10[1],wavelengthRangeInMUSELog10[2], by = 0.0001)
length(waveoutManual)
waveoutManual = seq(2, 9.35, by = 0.01)
testSED3 = ProSpectSED(massfunc = massfunc_dtau,
                       mSFR = 100,
                       mpeak = 10,
                       tau_birth = 1.0,
                       z = 0.0,
                       tau_screen = 0.3,
                       filtout=filtersData,
                       sparse=1,
                       speclib=BC03hr,
                       emission=FALSE
                       )
plot(testSED3, type="lum")
plot(testSED3$FinalFlux, log="xy", type="l")

plot(testSED3$FinalLum, log='xy', xlim=c(1e2,1e7), ylim=c(1e2,1e7),
     xlab='Wavelength (Ang)', ylab='Lum (Lsol/Ang)', type='l', lwd=5)
lines(testSED3$StarsUnAtten, col='blue', lty=2)
lines(testSED3$StarsAtten, col='green')
lines(testSED3$DustEmit, col='brown')
lines(testSED3$AGN, col='purple')
legend('topright',
       legend=c('Total Lum', 'Star Un-Atten', 'Stars Atten', 'Dust Emit', 'AGN'),
       col=c('black', 'blue', 'darkgreen', 'brown', 'purple'),
       lty=c(1,2,1,1,1),
       lwd=c(5,1,1,1,1)
)
# Probar cambiar el espectro de EMILES directamente para meterlo con un espectro reduciro (i.e. 3k-11k)

