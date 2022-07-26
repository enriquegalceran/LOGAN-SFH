setwd("~/Documents/GitHub/LOGAN-SFH")
library(RColorBrewer)
source("MUTANTS/LOGAN.R")
library(plyr)
library(ProSpect)
library(ggplot2)
library(plyr)         # splat(func)(c(var1, list_of_vars))
library(stringi)
library(dplyr)
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")
EMILESRecortado = readRDS(file="EMILESData/EMILESRecortado.rds")
'%!in%' <- function(x,y)!('%in%'(x,y))
EMILESCombined = readRDS(file="EMILESData/EMILESCombined.rds")





# 
# colvec = rev(rainbow(60))
# metal = 1
# if (metal==1) {maxage=50} else {maxage=60}
# for (age in 1:maxage){
#     if (age==1){
#         plot(EMILESCombined$Wave,
#              EMILESCombined$Zspec[[metal]][age, ]/max(EMILESCombined$Zspec[[1]][age, ])+(age-1),
#              type="l", log="x", col=colvec[age], ylim=c(0, maxage),
#              ylab="normalized flux", xlab="Wavelength",
#              main=paste0("Metallicity=", metal, " (", EMILESCombined$Z[[metal]], ")"))
#     } else {
#         lines(EMILESCombined$Wave,
#               EMILESCombined$Zspec[[metal]][age, ]/max(EMILESCombined$Zspec[[1]][age, ])+(age-1),
#               col=colvec[age])
#     }
# }
# abline(v=6500)
# abline(v=9000)
# for (h in 0:maxage){
#     abline(h=h, col="grey")
# }
# 
# plot(EMILESCombined$Wave,
#      EMILESCombined$Zspec[[1]][1,],
#      type="l", log="x", col=colvec[1],
#      ylab="flux", xlab="Wavelength",
#      main="Flux only for 6.3Myrs")
# abline(v=6500)
# abline(v=9000)


########################################3

filtersHST <- c("F275W", "F336W", "F438W", "F555W", "F814W")
filters <- list()
for (filter in filtersHST) {
    # TODO: see speclib regarding location of data files.
    filters[[filter]] <-
        read.table(
            paste0("FiltersHST/HST_WFC3_UVIS2.", filter, ".dat"),
            col.names = c("wave", "response")
        )
}

test_object = SFHfunc2(massfunc=massfunc_snorm_burst, mSFR=2.3911, mpeak=12.95,
        mperids=1.04, mskew=-0.298, mburstage=0.667, mburst=0.003,
        Z=Zfunc_massmap_box, Zstart=1e-3, Zfinal=1e-2,
        z=1e-4, veldisp=20000, yield=0.3, forcemass=3393729579, speclib=EMILESRecortado,
        filters=filters, emission=TRUE, emission_scale = "SFR", sparse=1
        )


plot(test_object$flux$wave, test_object$flux$flux, type="l", log="x")

a <- test_object$flux$wave[1:(length(test_object$flux$wave) - 1)] - test_object$flux$wave[2:length(test_object$flux$wave)]
plot(test_object$flux$wave[1:(length(test_object$flux$wave) - 1)], a, type="l", log="x")
plot(test_object$flux$wave[22:1057], a[22:1057], type="l", log="x")
plot(test_object$flux$wave[22:1057], test_object$flux$flux[22:1057], type="l", log="xy")


SFHfunc2 <- function (massfunc = massfunc_b5, forcemass = FALSE, agescale = 1, 
          stellpop = "BC03lr", speclib = NULL, tau_birth = 1, tau_screen = 0.3, 
          pow_birth = -0.7, pow_screen = -0.7, filters = "all", Z = 5, 
          emission = FALSE, veldisp = 50, emission_scale = "FUV", 
          escape_frac = 1 - emission, Ly_limit = 911.8, LKL10 = NULL, 
          z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, 
          ref, outtype = "mag", sparse = 5, intSFR = FALSE, unimax = 1.38e+10, 
          agemax = NULL, LumDist_Mpc = NULL, Eb = 0, L0 = 2175.8, 
          LFWHM = 470, ...) 
{
    dots = list(...)
    massfunc_args = dots[names(dots) %in% names(formals(massfunc))]
    if (stellpop == "BC03lr") {
        if (is.null(speclib)) {
            BC03lr = NULL
            data("BC03lr", envir = environment())
            speclib = BC03lr
        }
    }
    else if (stellpop == "BC03hr") {
        if (is.null(speclib)) {
            BC03hr = NULL
            data("BC03hr", envir = environment())
            speclib = BC03hr
        }
    }
    else if (stellpop == "EMILES") {
        if (is.null(speclib)) {
            EMILES = NULL
            data("EMILES", envir = environment())
            speclib = EMILES
        }
    }
    else if (stellpop == "BPASS") {
        if (is.null(speclib)) {
            BPASS = NULL
            data("BPASS", envir = environment())
            speclib = BPASS
        }
    }
    if (any(speclib$Age <= 1e+07)) {
        birthcloud = max(which(speclib$Age <= 1e+07))
    }
    else {
        birthcloud = 1
    }
    if (!is.null(filters)) {
        if (filters[1] == "all" | filters[1] == "GAMA") {
            filters = c("FUV_GALEX", "NUV_GALEX", "u_VST", "g_VST", 
                        "r_VST", "i_VST", "Z_VISTA", "Y_VISTA", "J_VISTA", 
                        "H_VISTA", "K_VISTA", "W1_WISE", "W2_WISE", 
                        "W3_WISE", "W4_WISE", "P100_Herschel", "P160_Herschel", 
                        "S250_Herschel", "S350_Herschel", "S500_Herschel")
        }
    }
    if (!is.function(Z)) {
        if (Z%%1 != 0) {
            tempZfunc = function(age, Z, ...) {
                rep(Z, length(age))
            }
            formals(tempZfunc)$Z = Z
            Z = tempZfunc
        }
    }
    agevec = speclib$Age * agescale
    if (is.function(Z)) {
        dots = list(...)
        Z_args = dots[names(dots) %in% names(formals(Z))]
        Zvec = do.call("Z", c(list(agevec), list(massfunc = massfunc), 
                              Z_args, massfunc_args))
        Zlist = interp_param(Zvec, speclib$Z, log = TRUE)
        Zwmat = matrix(0, length(speclib$Age), length(speclib$Z))
        Zwmat[cbind(1:length(speclib$Age), Zlist[, "ID_hi"])] = Zlist[, 
                                                                      "wt_hi"]
        Zwmat[cbind(1:length(speclib$Age), Zlist[, "ID_lo"])] = Zlist[, 
                                                                      "wt_lo"]
        Zuse = which(colSums(Zwmat) > 0)
        Zdoweight = TRUE
    }
    else {
        Zuse = Z
        Zvec = rep(speclib$Z[Zuse], length(speclib$Age))
        Zdoweight = FALSE
    }
    wave_lum = speclib$Wave
    if (sparse > 1) {
        sparse = seq(1, dim(speclib$Zspec[[1]])[2], by = sparse)
        for (i in Zuse) {
            speclib$Zspec[[i]] = speclib$Zspec[[i]][, sparse]
        }
        wave_lum = wave_lum[sparse]
    }
    if (intSFR) {
        massvec = rep(0, length(agevec))
        for (i in 1:length(agevec)) {
            tempint = try(do.call("integrate", c(list(f = massfunc, 
                                                      lower = speclib$AgeBins[i] * agescale, upper = speclib$AgeBins[i + 
                                                                                                                         1] * agescale), massfunc_args))$value, silent = TRUE)
            if (class(tempint) == "try-error") {
                massvec[i] = 0
            }
            else {
                massvec[i] = tempint
            }
        }
        if (sum(massvec, na.rm = TRUE) == 0) {
            massvec = rep(0, length(agevec))
            for (i in 1:length(agevec)) {
                tempint = try(do.call("integral", c(list(f = massfunc, 
                                                         xmin = speclib$AgeBins[i] * agescale, xmax = speclib$AgeBins[i + 
                                                                                                                          1] * agescale), massfunc_args)), silent = TRUE)
                if (class(tempint) == "try-error") {
                    massvec[i] = 0
                }
                else {
                    massvec[i] = tempint
                }
            }
        }
    }
    else {
        massvec = do.call("massfunc", c(list(agevec), massfunc_args)) * 
            speclib$AgeWeights
    }
    if (unimax != FALSE) {
        if (is.null(agemax)) {
            agemax = unimax - ProSpect::cosdistTravelTime(z = z, H0 = H0, 
                                                OmegaM = OmegaM, OmegaL = OmegaL, ref = ref) * 
                1e+09
        }
    }
    if (!is.null(agemax)) {
        massvec[speclib$Age > agemax] = 0
    }
    if (any(massvec < 0)) {
        stop("Supplied massfunc cannot create negative SFR!")
    }
    if (forcemass == FALSE) {
        forcescale = 1
        masstot = sum(massvec)
    }
    else {
        masstot = sum(massvec)
        forcescale = forcemass/masstot
        massvec = massvec * forcescale
        masstot = forcemass
    }
    if (length(Zuse) > 1) {
        lum = rep(0, length(wave_lum))
        for (Zid in Zuse) {
            lum = lum + colSums(speclib$Zspec[[Zid]] * massvec * 
                                    Zwmat[, Zid])
            if (any(escape_frac < 1)) {
                if (length(Ly_limit) == 1) {
                    speclib$Zspec[[Zid]][, wave_lum < Ly_limit] = speclib$Zspec[[Zid]][, 
                                                                                       wave_lum < Ly_limit] * escape_frac
                }
                else {
                    for (i in 1:(length(Ly_limit) - 1)) {
                        sel = which(wave_lum < Ly_limit[i] & wave_lum > 
                                        Ly_limit[i + 1])
                        speclib$Zspec[[Zid]][, sel] = speclib$Zspec[[Zid]][, 
                                                                           sel] * escape_frac[i]
                    }
                    sel = which(wave_lum < Ly_limit[length(Ly_limit)])
                    speclib$Zspec[[Zid]][, sel] = speclib$Zspec[[Zid]][, 
                                                                       sel] * escape_frac[length(Ly_limit)]
                }
            }
        }
    }
    else {
        lum = colSums(speclib$Zspec[[Zuse]] * massvec)
        if (any(escape_frac < 1)) {
            if (length(Ly_limit) == 1) {
                speclib$Zspec[[Zuse]][, wave_lum < Ly_limit] = speclib$Zspec[[Zuse]][, 
                                                                                     wave_lum < Ly_limit] * escape_frac
            }
            else {
                for (i in 1:(length(Ly_limit) - 1)) {
                    sel = which(wave_lum < Ly_limit[i] & wave_lum > 
                                    Ly_limit[i + 1])
                    speclib$Zspec[[Zuse]][, sel] = speclib$Zspec[[Zuse]][, 
                                                                         sel] * escape_frac[i]
                }
                sel = which(wave_lum < Ly_limit[length(Ly_limit)])
                speclib$Zspec[[Zuse]][, sel] = speclib$Zspec[[Zuse]][, 
                                                                     sel] * escape_frac[length(Ly_limit)]
            }
        }
    }
    lumtot_unatten = sum(.qdiff(wave_lum) * lum)
    lum_unatten = lum
    if (tau_birth != 0) {
        lum = rep(0, length(wave_lum))
        for (Zid in Zuse) {
            if (tau_birth != 0) {
                speclib$Zspec[[Zid]][1:birthcloud, ] = speclib$Zspec[[Zid]][1:birthcloud, 
                ] * rep(CF_birth(wave_lum, tau = tau_birth, 
                                 pow = pow_birth), each = birthcloud)
            }
            if (Zdoweight) {
                lum = lum + colSums(speclib$Zspec[[Zid]] * massvec * 
                                        Zwmat[, Zid])
            }
            else {
                lum = colSums(speclib$Zspec[[Zid]] * massvec)
            }
        }
        lumtot_birth = lumtot_unatten - sum(.qdiff(wave_lum) * 
                                                lum)
    }
    else {
        lumtot_birth = 0
    }
    if (emission) {
        if (emission_scale == "FUV") {
            if (length(Ly_limit) == 1 | all(escape_frac == escape_frac[1])) {
                sel = which(wave_lum < Ly_limit[1])
                All_lum = (1 - escape_frac) * sum(.qdiff(wave_lum[sel]) * 
                                                      lum_unatten[sel])
            }
            else {
                All_lum = 0
                for (i in 1:(length(Ly_limit) - 1)) {
                    sel = which(wave_lum < Ly_limit[i] & wave_lum > 
                                    Ly_limit[i + 1])
                    All_lum = All_lum + (1 - escape_frac[i]) * 
                        (Ly_limit[i] - Ly_limit[i + 1]) * mean(lum_unatten[sel])
                }
                sel = which(wave_lum < Ly_limit[length(Ly_limit)])
                All_lum = All_lum + (1 - escape_frac[length(Ly_limit)]) * 
                    sum(.qdiff(wave_lum[sel]) * lum_unatten[sel])
            }
            emission_input = list(All_lum = All_lum, veldisp = veldisp, 
                                  Z = Zvec[1])
            emissionadd_unatten = emissionLines(All_lum = All_lum, 
                                                veldisp = veldisp, Z = Zvec[1])
        }
        else if (emission_scale == "SFR") {
            SFRburst_emission = (1 - escape_frac) * do.call("integrate", 
                                                            c(list(f = massfunc, lower = 0, upper = 1e+07), 
                                                              massfunc_args))$value * forcescale/1e+07
            emission_input = list(SFR = SFRburst_emission, veldisp = veldisp, 
                                  Z = Zvec[1])
            emissionadd_unatten = emissionLines(SFR = SFRburst_emission, 
                                                veldisp = veldisp, Z = Zvec[1])
        }
        else {
            stop("emission_scale must be one of SFR or FUV!")
        }
        emissionadd_atten = emissionadd_unatten
        emissionadd_atten$lum = emissionadd_atten$lum * CF_birth(emissionadd_atten$wave, 
                                                                 tau = tau_birth, pow = pow_birth)
        lumtot_emission_unatten = sum(.qdiff(emissionadd_unatten$wave) * 
                                          emissionadd_unatten$lum)
        lumtot_emission_atten = sum(.qdiff(emissionadd_atten$wave) * 
                                        emissionadd_atten$lum)
        lumtot_birth = lumtot_birth - lumtot_emission_unatten + 
            lumtot_emission_atten
        if (lumtot_birth < 0) {
            lumtot_birth = 0
        }
        lum = addspec(wave_lum, lum, emissionadd_atten$wave, 
                      emissionadd_atten$lum)
        lum_unatten = approxfun(log10(wave_lum), lum_unatten)(log10(lum$wave))
        wave_lum = lum$wave
        lum = lum$flux
    }
    else {
        emission_input = NULL
    }
    if (tau_screen != 0) {
        lum = lum * CF_screen(wave_lum, tau = tau_screen, pow = pow_screen, 
                              Eb = Eb, L0 = L0, LFWHM = LFWHM)
        lumtot_screen = (lumtot_unatten - lumtot_birth) - sum(.qdiff(wave_lum) * 
                                                                  lum)
        if (lumtot_screen < 0) {
            lumtot_screen = 0
        }
    }
    else {
        lumtot_screen = 0
    }
    SFRburst = do.call("integrate", c(list(f = massfunc, lower = 0, 
                                           upper = 1e+08), massfunc_args))$value * forcescale/1e+08
    lumtot_atten = lumtot_unatten - sum(.qdiff(wave_lum) * lum)
    if (z < 0 | is.null(filters)) {
        out = NULL
        flux = NULL
    }
    else {
        if (z > 0) {
            flux = Lum2Flux(wave = wave_lum, lum = lum, z = z, 
                            H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref, 
                            LumDist_Mpc = LumDist_Mpc)
            if (!is.null(outtype)) {
                out = photom_flux(flux, outtype = outtype, filters = filters)
                if (is.list(filters)) {
                    cenout = {
                    }
                    for (i in filters) {
                        cenout = c(cenout, cenwavefunc(i))
                    }
                    if (all(!is.null(names(filters)))) {
                        out = data.frame(filter = names(filters), 
                                         cenwave = cenout, out = out)
                    }
                    else {
                        out = data.frame(filter = NA, cenwave = cenout, 
                                         out = out)
                    }
                }
                else {
                    cenwave = NULL
                    data("cenwave", envir = environment())
                    out = data.frame(cenwave[match(filters, cenwave$filter), 
                    ], out = out)
                }
            }
            else {
                out = NULL
            }
        }
        else {
            flux = cbind(wave = wave_lum, flux = lum * .lsol_to_absolute)
            if (!is.null(outtype)) {
                out = photom_flux(flux, outtype = outtype, filters = filters)
                if (is.list(filters)) {
                    cenout = {
                    }
                    for (i in filters) {
                        cenout = c(cenout, cenwavefunc(i))
                    }
                    if (all(!is.null(names(filters)))) {
                        out = data.frame(filter = names(filters), 
                                         cenwave = cenout, out = out)
                    }
                    else {
                        out = data.frame(filter = NA, cenwave = cenout, 
                                         out = out)
                    }
                }
                else {
                    cenwave = NULL
                    data("cenwave", envir = environment())
                    out = data.frame(cenwave[match(filters, cenwave$filter), 
                    ], out = out)
                }
            }
            else {
                out = NULL
            }
        }
    }
    SFR = massvec/speclib$AgeWeights
    if (!is.null(agemax)) {
        sel = which(speclib$Age <= agemax)
        agevec = agevec[sel]
        massvec = massvec[sel]
        SFR = SFR[sel]
        Zvec = Zvec[sel]
    }
    return(list(flux = flux, out = out, wave_lum = wave_lum, 
                lum_unatten = lum_unatten, lum_atten = lum, lumtot_unatten = lumtot_unatten, 
                lumtot_atten = lumtot_atten, lumtot_birth = lumtot_birth, 
                lumtot_screen = lumtot_screen, agevec = agevec, SFR = SFR, 
                masstot = masstot, massvec = massvec, M2L = masstot/lumtot_unatten, 
                SFRburst = SFRburst, Zvec = Zvec, emission_input = emission_input))
}
