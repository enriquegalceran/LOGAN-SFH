# Change pwd and load functions + libraries
# setwd("~/GitHub/Prospector")
setwd("~/Documents/GitHub/LOGAN-SFH")
source("FuncionesAuxiliaresYLibrerias.R")

data("BC03hr")
data("EMILES")


################################################################
#######################      FILTERS     #######################
################################################################

# Load the different HST filters into a
# HST filters: F275W (NUV), F336W (U), F438W (B), F555W (V), F814W (I)
filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filtersData <- list()
for (filter in filtersHST) {
  filtersData[[filter]] <-
    read.table(
      paste("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat", sep = ""),
      col.names = c("wave", "response")
    )
}

plotFilters(
  filtersHST,
  plottitle = "HST Filters",
  fromData = filtersData,
  sizeLegend = 1,
  setyLimTo1 = FALSE,
  singleLegend = FALSE,
  filternames = c("F275W(NUV)", "F336W(U)", "F438W(B)",
                  "F555W(V)", "F814W(I)")
)

# How to output TO A VARIABLE a single column from the filtersDatalist
# s <- "F275W"
# test <- filtersData[s]
# df <- data.frame(test)
# colnames(df) <- c("wave", "response")



################################################################
#######################        SFH       #######################
################################################################

# SFH for the moment will be set to dtau
# Remember that age = look-back time
starFormationTmax = 10
starFormationPeak = 200
SFHdtau <- function(x)
  massfunc_dtau(x,
                mSFR = starFormationPeak,
                mpeak = starFormationTmax,
                magemax = 10)
# Plot SFH
curve(SFHdtau, 0, 14e9, xlab = 'Age/yr', ylab = 'SFR / Msol/yr')



################################################################
#######################      BC03hr      #######################
################################################################

# Generate Dataset using BC03 (which we know is going to work)
SFHOutputBC03hr = SFHfunc(
  massfunc = massfunc_dtau,
  mSFR = starFormationPeak,
  mpeak = starFormationTmax,
  tau_birth = 1.0,
  z = 0.0,
  tau_screen = 0.3,
  sparse = 1,
  filters = filtersData,
  emission = TRUE,
  stellpop = "BC03hr"
)

plotonlyunatten(SFHOutputBC03hr, "BC03highres")
plotonlyunatten(SFHOutputBC03hr, "BC03highres", setxLim = "limit")



################################################################
#######################      BC03lr      #######################
################################################################

SFHOutputBC03lr = SFHfunc(
  massfunc = massfunc_dtau,
  mSFR = starFormationPeak,
  mpeak = starFormationTmax,
  tau_birth = 1.0,
  z = 0.0,
  tau_screen = 0.3,
  sparse = 1,
  filters = filtersData,
  emission = TRUE,
  magemax = 13.8,
  stellpop = "BC03lr"
)

plotonlyunatten(SFHOutputBC03lr, "BC03lowres")
plot(SFHOutputBC03lr$agevec, SFHOutputBC03lr$SFR, type = "l")


################################################################
#######################      EMILES      #######################
################################################################


SFHOutputEMILES = SFHfunc(
  massfunc = massfunc_dtau,
  mSFR = starFormationPeak,
  mpeak = starFormationTmax,
  magemax = 13.8,
  tau_birth = 1.0,
  z = 0.0,
  tau_screen = 0.3,
  sparse = 1,
  filters = filtersData,
  emission = TRUE,
  stellpop = "EMILES",
  emission_scale = 'SFR',
  Z = Zfunc_massmap_box,
  yield = 0.035,
  Zfinal = 0.05
)
SFHOutputEMILESNoEm = SFHfunc(
  massfunc = massfunc_dtau,
  mSFR = starFormationPeak,
  mpeak = starFormationTmax,
  magemax = 13.8,
  tau_birth = 1.0,
  z = 0.0,
  tau_screen = 0.3,
  sparse = 1,
  filters = filtersData,
  emission = FALSE,
  stellpop = "EMILES",
  emission_scale = 'SFR',
  Z = Zfunc_massmap_box,
  yield = 0.035,
  Zfinal = 0.05
)



plotonlyunatten(SFHOutputEMILES, "EMILES")

plotonlyflux(SFHOutputEMILES, "EMILES")



###################
# function(x1,x2,x3,x4, ...){
#   dots = list(...)
# }

exportFITSfile(SFHOutputEMILES, "EMILESOutput")

readDataEMILESParamt <- readFITSfile("Paramt_EMILESOutput")
readDataEMILESSpectr <- readFITSfile("Spectr_EMILESOutput")
readDataEMILESFilter <- readFITSfile("Filter_EMILESOutput")



ProSpectSED()



