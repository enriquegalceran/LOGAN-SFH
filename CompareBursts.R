# Compare spectra for bursts

# Initialize variables
.pardefault <- par()
setwd("~/Documents/GitHub/LOGAN-SFH")
library(ProSpect)
library(RColorBrewer)
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")

# filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
# filters <- list()
# for (filter in filtersHST) {
#   filters[[filter]] <-
#     read.table(
#       paste("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat", sep = ""),
#       col.names = c("wave", "response")
#     )
# }



age = EMILESCombined$Age



massParams3a = list(
  snorm_burst=list(
    name="snorm_burst",
    func=massfunc_snorm_burst,
    mskew=seq(-0.25, 0.25, 0.25),
    mperiod=seq(0.25, 1, 0.25),
    mpeak=seq(9, 12, 1),
    mburstage=seq(0.01, 0.09, 0.02),
    mburst=seq(1, 5, 1)
  )
)

func = massfunc_snorm_burst
common.parameters = list()

y <- func(age, mburst=4, mburstage = 0.1)
plot(age, y, type="l", log="x")

massfunc_cut <- function(age, mcut, massfunc_original, ...){
  dots = list(...)
  massfunc_args = dots[names(dots) %in% names(formals(massfunc_original))]
  
  mass_out = do.call("massfunc_original", c(list(age), massfunc_args))
  
  for (ageIdx in 1:length(age)){
    if (age[ageIdx] > mcut * 1e+09){
      mass_out[ageIdx] = 0
    }
  }
  return(mass_out)
}


y2 <- massfunc_cut(age, mcut=0.1, massfunc_snorm_burst, mburst=4, mburstage = 0.1)
plot(age, y2, type="l", log="x")


SFHfunc()

names(formals(massfunc_cut))

formals(massfunc_cut) <- c(formals(massfunc_cut), list(mskew=0.1))
formals(massfunc_cut)








