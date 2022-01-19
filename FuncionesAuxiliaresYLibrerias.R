# Librerías a cargar
library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)
library(latex2exp)
library(FITSio)
library(plotrix)


# Data Output
library("rio")
library("readr")
library("data.table")
library("feather")

source('Constants.R')
.pardefault <- par()

# Obtener t_0 y mu para las funciones log normal
getMuFromT0 <- function(t0, tau)
  log(t0) + 1 / tau
getT0FromMu <- function(mu, tau)
  exp(mu - 1 / t)

# Wavelength Range in MUSE
wavelengthRangeInMUSE <- c(4750, 9350)  # Ang

findClosestInVector <- function(number, vector, returnIndex = FALSE) {
  if (returnIndex) {
    return(which.min(abs(vector - number)))
  } else {
    return(vector[which.min(abs(vector - number))])
  }
}

# Conseguir memoria ram de las variables y objetos
sort(sapply(ls(), function(x) {
  object.size(get(x))
}))

wavelength2rgb <- function(l) {
  r = 0
  g = 0
  b = 0
  if ((l >= 400.0) && (l < 410.0)) {
    t = (l - 400.0) / (410.0 - 400.0)
    r =    +(0.33 * t) - (0.20 * t * t)
  } else if ((l >= 410.0) && (l < 475.0)) {
    t = (l - 410.0) / (475.0 - 410.0)
    r = 0.14         - (0.13 * t * t)
  } else if ((l >= 545.0) && (l < 595.0)) {
    t = (l - 545.0) / (595.0 - 545.0)
    r =    +(1.98 * t) - (t * t)
  } else if ((l >= 595.0) && (l < 650.0)) {
    t = (l - 595.0) / (650.0 - 595.0)
    r = 0.98 + (0.06 * t) - (0.40 * t * t)
  } else if ((l >= 650.0) && (l < 700.0)) {
    t = (l - 650.0) / (700.0 - 650.0)
    r = 0.65 - (0.84 * t) + (0.20 * t * t)
  }
  if ((l >= 415.0) && (l < 475.0)) {
    t = (l - 415.0) / (475.0 - 415.0)
    g =             +(0.80 * t * t)
  } else if ((l >= 475.0) && (l < 590.0)) {
    t = (l - 475.0) / (590.0 - 475.0)
    g = 0.8 + (0.76 * t) - (0.80 * t * t)
  } else if ((l >= 585.0) && (l < 639.0)) {
    t = (l - 585.0) / (639.0 - 585.0)
    g = 0.84 - (0.84 * t)
  }
  if ((l >= 400.0) && (l < 475.0)) {
    t = (l - 400.0) / (475.0 - 400.0)
    b =    +(2.20 * t) - (1.50 * t * t)
  } else if ((l >= 475.0) && (l < 560.0)) {
    t = (l - 475.0) / (560.0 - 475.0)
    b = 0.7 - (t) + (0.30 * t * t)
  }
  return(c(r, g, b))
}


# Ayuda para dibujar gráficas
plotmanydata <- function(df,
                         name = NULL,
                         setxLim = "lines",
                         waveRangeMUSE = wavelengthRangeInMUSE) {
  if (setxLim == "limit")
    xLim = waveRangeMUSE
  else
    xLim = ""
  
  plot(
    df$flux[, 1],
    df$flux[, 2],
    type = "l",
    log = "xy",
    main = paste("Flux for ", name),
    xlab = "Wavelength (\uc5)",
    ylab = "Flux (erg/s/cm^2/\uc5)",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
  
  plot(
    df$out$cenwave,
    df$out$out,
    log = "xy",
    main = paste("Requested Output for central wavelength (", name, ")"),
    xlab = "Wavelength (\uc5)",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
  
  plot(
    df$wave_lum,
    df$lum_unatten,
    type = "l",
    log = "xy",
    main = paste("Unattenuated Stellar Luminosity (", name, ")"),
    xlab = "Wavelength (\uc5)",
    ylab = "L_sol",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
  
  plot(
    df$wave_lum,
    df$lum_atten,
    type = "l",
    log = "xy",
    main = paste("Attenuated Stellar Luminosity (", name, ")"),
    xlab = "Wavelength (\uc5)",
    ylab = "L_sol",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
  
  plot(
    df$agevec,
    df$SFR,
    main = paste("Star Formation Rate (", name, ")"),
    xlab = "ages",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
  
  plot(
    df$agevec,
    df$massvec,
    main = paste("total mass formed in each age window (", name, ")"),
    xlab = "ages",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
  
  plot(
    df$agevec,
    df$Zvec,
    main = paste("Metallicity(?) (", name, ")"),
    xlab = "ages",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
}


plotonlyflux <- function(df,
                         name = NULL,
                         setxLim = "lines",
                         waveRangeMUSE = wavelengthRangeInMUSE) {
  if (setxLim == "limit")
    xLim = waveRangeMUSE
  else
    xLim = NULL
  
  plot(
    df$flux[, 1],
    df$flux[, 2],
    type = "l",
    log = "xy",
    main = paste("Flux for", name),
    xlab = "Wavelength (\uc5)",
    ylab = "Flux (erg/s/cm^2/\uc5)",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
}


plotonlyunatten <- function(df,
                            name = NULL,
                            setxLim = "lines",
                            waveRangeMUSE = wavelengthRangeInMUSE) {
  if (setxLim == "limit")
    xLim = waveRangeMUSE
  else
    xLim = NULL
  
  plot(
    df$wave_lum,
    df$lum_unatten,
    type = "l",
    log = "xy",
    main = paste("Unattenuated Stellar Luminosity (", name, ")"),
    xlab = "Wavelength (\uc5)",
    ylab = "L_sol",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
}


plotonlyatten <- function(df,
                          name = NULL,
                          setxLim = "lines",
                          waveRangeMUSE = wavelengthRangeInMUSE) {
  if (setxLim == "limit")
    xLim = waveRangeMUSE
  else
    xLim = NULL
  
  plot(
    df$wave_lum,
    df$lum_atten,
    type = "l",
    log = "xy",
    main = paste("Attenuated Stellar Luminosity (", name, ")"),
    xlab = "Wavelength (\uc5)",
    ylab = "L_sol",
    pch = 16,
    cex = 0.1,
    xlim = xLim
  )
  if (setxLim == "lines") {
    abline(v = wavelengthRangeInMUSE[1])
    abline(v = wavelengthRangeInMUSE[2])
  }
}


plotFilters <- function(filterlist,
                        plottitle = NULL,
                        margin = 0.1,
                        log = "",
                        singleLegend = TRUE,
                        sizeLegend = 0.5,
                        wavelengthcolor = FALSE,
                        fromData = NULL,
                        setyLimTo1 = TRUE,
                        setxLim = "lines",
                        legendPlacement = "topright",
                        waveRangeMUSE = wavelengthRangeInMUSE,
                        filternames = NULL) {
  # Get minimum and maximum value of the filters wavelengths
  minTmp = Inf
  maxTmp = -Inf
  maxWavelength <- c()
  maxResponse <- c()
  for (filter in filterlist) {
    if (is.null(fromData)) {
      tmp <- get(paste('filt_', filter, sep = ""))
    } else {
      tmp <- fromData[filter]
      tmp <- data.frame(tmp)
      colnames(tmp) <- c("wave", "response")
    }
    
    if (min(tmp$wave) < minTmp) {
      minTmp = min(tmp$wave)
    }
    if (max(tmp$wave) > maxTmp) {
      maxTmp = max(tmp$wave)
    }
    
    # Calculate position of maximum response in filter
    if (length(tmp$wave[tmp$response == max(tmp$response)]) == 1) {
      maxWavelength <-
        c(maxWavelength, tmp$wave[tmp$response == max(tmp$response)])
      maxResponse <- c(maxResponse, max(tmp$response))
    } else {
      tmpwave = tmp$wave[tmp$response == max(tmp$response)]
      maxWavelength <- c(maxWavelength, tmpwave[1])
      maxResponse <-
        c(maxResponse, tmp$response[tmp$wave == tmpwave[1]])
    }
  }
  
  if (setyLimTo1)
    yLim = c(0, 1 + 0.02 * !singleLegend)
  else
    yLim = c(0, max(maxResponse) + 0.02 * !singleLegend)
  
  # setxLim valid input: lines, limit, none/NULL
  if (setxLim == "limit")
    xLim = waveRangeMUSE
  else
    xLim = c(minTmp * (1 - margin), maxTmp * (1 + margin))
  
  # If filtername != NULL, assert length
  if (!is.null(filternames)) {
    stopifnot(length(filternames) == length(filterlist))
  } else {
    filternames <- filterlist
  }
  
  
  # Colors
  if (wavelengthcolor == FALSE || wavelengthcolor == TRUE) {
    # Sort color according to the wavelength of the maximum response.
    # TODO: Can be further improved by obtaining the middle of the integrated value.
    colvec = rev(rainbow(length(filterlist), end = 2 / 3))
    colvec = colvec[order(maxWavelength)]
  }
  # else {
  # In development
  # colvec = c()
  # for (i in 1:length(filterlist)){
  #   print("-------")
  #   print(i)
  #   print(maxWavelength[i])
  #   waveRGB = wavelength2rgb(maxWavelength[i]/10)
  #   print(waveRGB)
  #   tmp = rgb(waveRGB[1], waveRGB[2], waveRGB[3], maxColorValue = 1)
  #   print(tmp)
  # }
  #
  # }
  
  # Plotting
  for (i in 1:length(filterlist)) {
    if (is.null(fromData)) {
      tmpObject = get(paste('filt_', filterlist[i], sep = ""))
    } else {
      tmpObject <- fromData[filterlist[i]]
      tmpObject <- data.frame(tmpObject)
      colnames(tmpObject) <- c("wave", "response")
    }
    
    if (i == 1) {
      plot(
        tmpObject,
        type = 'l',
        col = colvec[1],
        xlim = xLim,
        ylim = yLim,
        log = log,
        xlab = 'Wave / Ang (\uc5)',
        ylab = 'Response',
        main = plottitle
      )
    } else {
      lines(tmpObject, type = 'l', col = colvec[i])
    }
    if (singleLegend == FALSE) {
      text(maxWavelength[i], maxResponse[i] + 0.015, filternames[i])
    }
  }
  
  # If setxLim=="lines", draw vertical lines representing the range of MUSE
  if (setxLim == "lines") {
    abline(v = waveRangeMUSE[1], col = "black")
    abline(v = waveRangeMUSE[2], col = "black")
  }
  
  
  # Plot general
  if (singleLegend == TRUE)
    legend(
      legendPlacement,
      legend = filternames,
      col = colvec,
      lwd = 2,
      cex = sizeLegend
    )
  
}


testplotHST <- function() {
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
    fromData = filtersData,
    sizeLegend = 1,
    setyLimTo1 = FALSE
  )
}


exportData <- function(df,
                       filename,
                       foldername = "OutputFolder",
                       fileprefix = "",
                       fileformat = ".csv") {
  spectraOutput = data.frame(df$flux)
  filtersOutput = data.frame(df$out)
  parametOutput = data.frame(
    agevec = df$agevec,
    massvec = df$massvec,
    SFR = df$SFR,
    Z = df$Zvec
  )
  
  filedirectory = file.path(getwd(), foldername)
  # If output folder does not exist, generate it
  dir.create(filedirectory, showWarnings = FALSE)
  
  filename = paste(fileprefix, filename, fileformat, sep = "")
  
  cat("Exporting files in ", filedirectory, "/\n", sep = "")
  cat("Exporting Spectr_", filename, " ...\n", sep = "")
  export(spectraOutput,
         paste(filedirectory, "/Spectr_", filename, sep = ""))
  cat("Exporting Filter_", filename, " ...\n", sep = "")
  export(filtersOutput,
         paste(filedirectory, "/Filter_", filename, sep = ""))
  cat("Exporting Paramt_", filename, " ...\n", sep = "")
  export(parametOutput,
         paste(filedirectory, "/Paramt_", filename, sep = ""))
}


exportFITSfile <- function(df,
                           filename,
                           foldername = "OutputFolder",
                           fileprefix = "") {
  addFiltersToHeader <- function(header, filterlist) {
    header <- addComment("List of filters analysed.", header = header)
    for (i in 1:length(filterlist)) {
      header <-
        addKwv(paste("filter", i, sep = ""), filterlist[i], header = header)
    }
    return(header)
  }
  
  addColumnInformation <-
    function(header, columnnames, columnform = NULL) {
      header <- addComment("Columns information:", header = header)
      for (i in 1:length(columnnames)) {
        header <-
          addKwv(paste("TTYPE", i, sep = ""), columnnames[i], header = header)
        if (!is.null(columnform))
          header <-
            addKwv(paste("TFORM", i, sep = ""), columnform[i], header = header)
      }
      return(header)
    }
  
  spectraOutput = data.frame(df$flux)
  filtersOutput = data.frame(df$out)
  parametOutput = data.frame(
    agevec = df$agevec,
    massvec = df$massvec,
    SFR = df$SFR,
    Z = df$Zvec
  )
  
  filedirectory = file.path(getwd(), foldername)
  # If output folder does not exist, generate it
  dir.create(filedirectory, showWarnings = FALSE)
  filename = paste(fileprefix, filename, ".fits", sep = "")
  
  
  cat("Exporting files in ", filedirectory, "/\n", sep = "")
  ############################################################
  ########## Flux
  matrixFlux = data.matrix(spectraOutput[1:2])
  
  # make main header
  header <- makeFITSimHdr(
    matrixFlux,
    c1 = "FITS file Flux",
    ctypen = c("wave", "flux"),
    cunitn = c("A", "erg/s/cm**2/A")
  )
  header <-
    addColumnInformation(header, attributes(matrixFlux)$dimnames[[2]])
  # Check if file exists and remove if necessary
  if (file.exists(paste(filedirectory, "/Spectr_", filename, sep = ""))) {
    #Delete file if it exists
    file.remove(paste(filedirectory, "/Spectr_", filename, sep = ""))
  }
  # Write file
  cat("Exporting Spectr_", filename, " ...\n", sep = "")
  writeFITSim(
    matrixFlux,
    file = paste(filedirectory, "/Spectr_", filename, sep = ""),
    header = header
  )
  
  
  ############################################################
  ########## Filters
  
  matrixFilters <- data.matrix(filtersOutput[, 2:3])
  
  header <- makeFITSimHdr(
    matrixFilters,
    c1 = "FITS file Filters",
    ctypen = c("cenwave", "out"),
    cunitn = c("A", "mag")
  )
  header <- addFiltersToHeader(header, df$out$filter)
  header <-
    addColumnInformation(header, attributes(matrixFilters)$dimnames[[2]])
  
  # Check if file exists and remove if necessary
  if (file.exists(paste(filedirectory, "/Filter_", filename, sep = ""))) {
    #Delete file if it exists
    file.remove(paste(filedirectory, "/Filter_", filename, sep = ""))
  }
  # Write file
  cat("Exporting Filter_", filename, " ...\n", sep = "")
  writeFITSim(
    matrixFilters,
    file = paste(filedirectory, "/Filter_", filename, sep = ""),
    header = header
  )
  
  ############################################################
  ########## Parameters
  
  matrixParamt <- data.matrix(parametOutput)
  
  header <- makeFITSimHdr(matrixParamt)
  header <-
    addColumnInformation(header, attributes(matrixParamt)$dimnames[[2]])
  # Check if file exists and remove if necessary
  if (file.exists(paste(filedirectory, "/Paramt_", filename, sep = "")))
    file.remove(paste(filedirectory, "/Paramt_", filename, sep = ""))
  
  # Write file
  cat("Exporting Paramt_", filename, " ...\n", sep = "")
  writeFITSim(
    matrixParamt,
    file = paste(filedirectory, "/Paramt_", filename, sep = ""),
    header = header
  )
}





readFITSfile <- function(filename,
                         foldername = "OutputFolder",
                         fileprefix = "") {
  # Update filename
  filedirectory = file.path(getwd(), foldername)
  filename = paste(filedirectory, "/", fileprefix, filename, sep = "")
  if (substr(filename, nchar(filename) - 5, nchar(filename)) != ".fits")
    filename <- paste(filename, ".fits", sep = "")
  
  # Read File
  dfRead <-  readFITS(filename)
  
  # Generate list of column names
  cNames <- c()
  for (i in 1:length(dfRead$hdr))
    if (grepl("TTYPE", dfRead$hdr[i], fixed = TRUE))
      cNames <- c(cNames, dfRead$hdr[i + 1])
  
  # Read Subgroup
  dfData <- data.frame(dfRead$imDat)
  colnames(dfData) <- cNames
  
  return(dfData)
}




customWaveout <- function(fluxObject, waveout) {
  outObject = fluxObject
  
  waveoutLog = log10(waveout)
  waveLog = log10(fluxObject$wave_lum)
  unattenLog = log10(fluxObject$lum_unatten)
  attenLog = log10(fluxObject$lum_atten)
  
  outObject$wave_lum = waveout
  outObject$lum_unatten = 10 ^ approxfun(waveLog, unattenLog, rule = 2)(waveoutLog)
  outObject$lum_atten = 10 ^ approxfun(waveLog, attenLog, rule = 2)(waveoutLog)
  return(outObject)
}


interpolateToWaveout <- function(x1,
                                 y1,
                                 waveout,
                                 returnList = FALSE,
                                 offset = 0.5,
                                 n.points.integrate = 50,
                                 interpolate = FALSE,
                                 method.to.evaluate = "mean") {
  if (interpolate) {
    # If interpolate is TRUE, simply interpolate for the xnew (waveout).
    # This is NOT ADVISED if the NEW RESOLUTION IS SMALLER.
    waveoutL = log10(waveout)
    x1L = log10(x1)
    y1L = log10(y1)
    spect = 10 ^ approxfun(x1L, y1L, rule = 2)(waveoutL)
    if (!returnList) {
      return(spect)
    } else {
      return(list(wave = waveout, spect = spect))
    }
  } else {
    # If interpolate is FALSE, calculate the average values for the points for the new x
    
    # offset calculates where the bin limits will be placed:
    # the centerpoint c_i between two consecutive points (x_i and x_i+1, with x_i+1 > x_i) will be placed at
    # c_i <- (x_i+1 - x_i+1) * offset + x_i
    # the new evaluation will be:
    # y_i <- index in data where value corresponds to x=c_i
    # f_new(x_i) <- mean(data[y_i-1 : y_i])
    
    if (offset < 0 || offset > 1) {
      stop("'offset' needs to be within [0,1]")
    }
    
    # Initialize Y
    last.w = length(waveout)
    newy = numeric(last.w)
    
    # Iterate over every new waveout and define limits for integration
    for (i in 1:last.w) {
      # Separate between first, last and rest
      if (i == 1) {
        left  = waveout[1] -
                    (waveout[2] - waveout[1]) * (1 - offset)
        right = waveout[1] +
                    (waveout[2] - waveout[1]) * offset
        
      } else if (i == last.w) {
        left  = waveout[last.w - 1] + 
                    (waveout[last.w] - waveout[last.w - 1]) * offset
        right = waveout[last.w] + 
                    (waveout[last.w] - waveout[last.w -1]) * offset
        
      } else {
        left  = waveout[i - 1] +
                    (waveout[i] - waveout[i - 1]) * offset
        right = waveout[i]   +
                    (waveout[i + 1] - waveout[i]) * offset
      }
      
      # Once we know the limits between which we need to integrate,
      # interpolate a higher resolution (n.points.integrate) x vector
      # and interpolate the new values (using logs).
      tmp = interpolateToWaveout(x1,
                                 y1,
                                 seq(left, right, length.out = n.points.integrate),
                                 interpolate = TRUE)
      
      # Obtain the final value from this vector (mean/median)
      if (method.to.evaluate == "mean") {
        newy[i] = mean(tmp)
      } else if (method.to.evaluate == "median") {
        newy[i] = median(tmp)
      }
    }
    return(newy)
  }
}


qdiffCustom <- function (vec, pad0 = TRUE) {
  if (pad0) {
    return(c(0, vec[2:length(vec)] - vec[1:(length(vec) - 1)]))
  }
  else {
    return(vec[2:length(vec)] - vec[1:(length(vec) - 1)])
  }
}


pltComparison <- function(object1, object2, mode = NULL) {
  if (is.null(mode)) {
    cat("Mode has the following accepted values:\ndifference1\ndifference2")
  }
  print("inside Function")
  margin = 0.05
  # Calculate data and margins
  {
    limitsSpecLib = c(1680, 49999)
    insideRangeIdx1 = object1$wave_lum >= limitsSpecLib[1] &
      object1$wave_lum <= limitsSpecLib[2]
    insideRangeIdx2 = object2$wave_lum >= limitsSpecLib[1] &
      object2$wave_lum <= limitsSpecLib[2]
    
    x1 = object1$wave_lum[insideRangeIdx1]
    x2 = object2$wave_lum[insideRangeIdx2]
    y11 = object1$lum_unatten[insideRangeIdx1]
    y12 = object1$lum_atten[insideRangeIdx1]
    y21 = object2$lum_unatten[insideRangeIdx2]
    y22 = object2$lum_atten[insideRangeIdx2]
    
    # Margins
    xLim = c(min(c(x1, x2), na.rm = TRUE), max(c(x1, x2), na.rm = TRUE))
    yLim = c(min(c(y11, y12, y21, y22), na.rm = TRUE),
             max(c(y11, y12, y21, y22), na.rm = TRUE))
    xLim = log10(xLim)
    yLim = log10(yLim)
    xLim = xLim + c(-1,+1) * margin * (xLim[2] - xLim[1]) / 2
    yLim = yLim + c(-1,+1) * margin * (yLim[2] - yLim[1]) / 2
    xLim = 10 ^ xLim
    yLim = 10 ^ yLim
  }
  
  par(oma = c(1, 1, 4, 3))    # Outer margins
  par(mai = c(0, 0.5, 0, 0))    # internal margins
  layout(matrix(
    c(1, 1, 1, 1, 1, 2, 2),
    nrow = 7,
    ncol = 1,
    byrow = TRUE
  ))
  
  plot(
    object1$wave_lum,
    object1$lum_unatten,
    log = "xy",
    xlim = xLim,
    ylim = yLim,
    ylab = "Luminosity (Lsol)",
    type = "l",
    col = "blue",
    xaxt = "n"
  )
  lines(object1$wave_lum, object1$lum_atten, col = 'green')
  lines(object2$wave_lum, object2$lum_unatten, col = 'red')
  lines(object2$wave_lum, object2$lum_atten, col = 'maroon')
  abline(v = wavelengthRangeInMUSE[1], col = "black")
  abline(v = wavelengthRangeInMUSE[2], col = "black")
  legend(
    'topright',
    legend = c("Unatt-O", "Atten-O", "Unatt-Y+O", "Atten-Y+O", "MUSE"),
    col = c('blue', 'green', 'red', 'maroon', 'black'),
    lty = c(1, 1, 1, 1),
    cex = 0.75
  )
  
  par(mai = c(0.5, 0.5, 0, 0))
  if (mode == "difference1") {
    y21b = interpolateToWaveout(x2, y21, x1, returnList = TRUE)
    y22b = interpolateToWaveout(x2, y22, x1, returnList = TRUE)
    differenceUnatten = y11 - y21b
    differenceUnattenAbs = abs(differenceUnatten)
    reldifferenceUnatten = differenceUnattenAbs / y11 * 100
    signUnatten = differenceUnatten == differenceUnattenAbs
    colorDifference = 1 * signUnatten
    sectors = c(1)
    tmpSign = signUnatten[1]
    if (tmpSign) {
      colSectors = c("green")
    } else {
      colSectors = c("red")
    }
    for (i in 1:length(signUnatten)) {
      if (signUnatten[i] != tmpSign) {
        sectors = c(sectors, i)
        tmpSign = signUnatten[i]
        if (tmpSign) {
          colSectors = c(colSectors, "green")
        } else {
          colSectors = c(colSectors, "red")
        }
      }
    }
    sectors = c(sectors, i)
    
    plot(
      x1[1],
      reldifferenceUnatten[1],
      log = "xy",
      xlim = xLim,
      ylim = c(min(reldifferenceUnatten), max(reldifferenceUnatten)),
      xlab = "Wavelength (Ang)",
      ylab = "Relative (%)",
      type = "l",
      col = "blue"
    )
    legend(
      "topright",
      legend = c("positive", "negative"),
      lty = c(1, 1),
      col = c("green", "red")
    )
    
    for (i in 1:(length(sectors) - 1)) {
      lines(x1[sectors[i]:sectors[i + 1]],
            reldifferenceUnatten[sectors[i]:sectors[i + 1]],
            col = colSectors[i])
    }
  } else if (mode == "difference2") {
    y21b = interpolateToWaveout(x2, y21, x1, returnList = TRUE)
    y22b = interpolateToWaveout(x2, y22, x1, returnList = TRUE)
    differenceUnatten = y11 - y21b
    reldifferenceUnatten = differenceUnatten / y11 * 100
    differenceAtten = y11 - y21b
    reldifferenceAtten = differenceAtten / y11 * 100
    
    plot(
      x1,
      reldifferenceUnatten,
      log = "x",
      xlim = xLim,
      ylim = c(min(reldifferenceUnatten), max(reldifferenceUnatten)),
      xlab = "Wavelength (Ang)",
      ylab = "Relative difference (%)",
      type = "l",
      col = "blue"
    )
    lines(x1, reldifferenceAtten, col = "blue")
    abline(h = 0, col = "black")
  } else if (mode == "division") {
    y21b = interpolateToWaveout(x2, y21, x1, returnList = TRUE)
    y22b = interpolateToWaveout(x2, y22, x1, returnList = TRUE)
    differenceUnatten = log10(y21b / y11)
    differenceAtten = log10(y22b / y12)
    plot(
      x1,
      differenceUnatten,
      log = "x",
      xlim = xLim,
      ylim = c(min(
        c(differenceUnatten, differenceAtten)
      ), max(
        c(differenceUnatten, differenceAtten)
      )),
      xlab = "Wavelength (Ang)",
      ylab = "[(Y+O)/O]",
      type = "l",
      col = "blue"
    )
    lines(x1, differenceAtten, col = "green")
    abline(h = 0, col = "black")
    abline(v = wavelengthRangeInMUSE[1], col = "black")
    abline(v = wavelengthRangeInMUSE[2], col = "black")
    legend(
      "topright",
      legend = c("Unatten", "Atten"),
      col = c("blue", "green"),
      lty = c(1, 1)
    )
  } else {
    cat("Esta opción todavía no está considerada.")
  }
  
  par(.pardefault)
}


convertAny2decimal <- function(x, old = 16) {
  if (old > 62) {
    stop("Maximum order is 62.")
  } else if (old <= 0) {
    stop("Minimum order is 1.")
  }
  x <- toString(x)
  dic <- c(0:9, letters[1:26], LETTERS[1:26])
  tmp = 0
  tmpi = 0
  while (x != "") {
    t <- match(substr(x, nchar(x), nchar(x)), dic) - 1
    tmp = tmp + old ** tmpi * t
    tmpi = tmpi + 1
    x = substr(x, 1, nchar(x) - 1)
  }
  return(tmp)
}


convertDecimal2Any <- function(x, new = 16) {
  if (new > 62) {
    stop("Maximum order is 62.")
  } else if (new <= 0) {
    stop("Minimum order is 1.")
  }
  dic <- c(0:9, letters[1:26], LETTERS[1:26])
  tmp = ""
  while (x > 0) {
    tmp = paste0(dic[x %% new + 1], tmp)
    x = (x - x %% new) / new
  }
  return(tmp)
}


convertAny2Any <- function(x, old = 10, new = 16) {
  if (old > 62 | new > 62) {
    stop("Maximum order is 62.")
  } else if (old <= 0 | new <= 0) {
    stop("Minimum order is 1.")
  }
  dic <- c(0:9, letters[1:26], LETTERS[1:26])
  convertDecimal2Any(convertAny2decimal(x, old), new)
}


convertAgevecToOutputScale <- function(agevector,
                                       datavector,
                                       new_scale=NULL,
                                       max.age=13.8e9,
                                       n.splits=10,
                                       method="mean",
                                       return_scale=FALSE) {
  # Converts Agevector (based on the Library [EMILESCombined]) to a smaller vector (for the CNN output).
  if (new_scale == "defaultlog1" || is.null(new_scale)) {
    new_scale = c(0, 6.3, 7.9, 10, 12.6, 15.8, 20, 25.1, 31.6, 39.8, 50.1, 63.1, 70.8,
                  130, 255, 510, 900, 1800, 3200, 6400, 12600)
  } else if (new_scale == "defaultlog2") {
    new_scale = c(0, 10, 32, 50.5, 130, 255, 510, 900, 1800, 3200, 6400, 12600)
  } else if (new_scale == "lovell"){
    new_scale = c(0, 32, 68, 147, 316, 681, 1470, 3160, 12460)
  } else if (new_scale == "newlog") {
    # Separate into n.splits bins (log)
    agevec.below.max = agevector[agevector<max.age]
    log_agevec = log10(agevec.below.max)
    
    # Split homogeneously among the values
    age_split <- split(agevec.below.max, cut(seq_along(log_agevec), n.splits, labels = FALSE))
    
    # Keep the last value (this gives the max value for the bin)
    new_scale = c(0)
    for (i in 1:length(age_split)){
      new_scale = c(new_scale, age_split[[i]][length(age_split[[i]])]/1000000)
    }
  } # ToDo: Maybe a 3rd default value that uses 
  
  # Identify which values of the agevector should go in each bin. NA if a value is above the max value (las element in new_scale)
  separations <- cut(agevector, new_scale * 1000000, labels = FALSE, include.lowest = TRUE)
  
  # Separate the data according to the x value
  data_groups <- split(datavector, separations)
  age_groups <- split(agevector, separations)
  
  # Consolidate the data using mean or median
  out <- c()
  age <- c()
  width <- c()
  for (i in 1:length(data_groups)) {
    # First, get the data into the new bins
    if (method == "mean") {
      tmp <- mean(data_groups[[i]])
    } else if (method == "median") {
      tmp <- median(data_groups[[i]])
    }
    out <- c(out, tmp)
    
    width <- c(width, length(age_groups[[i]]))
    
    # Second, get the ages of the bins
    if (width[i] == 1){
      age <- c(age, age_groups[[i]])
    } else {
      # age <- c(age, mean(age_groups[[i]]))
      age <- c(age, 10^mean(log10(age_groups[[i]])))
    }
  }
  if (return_scale){
    output = list(age=age, data=out, width=width, new_scale=new_scale)
  } else {
    output = list(age=age, data=out, width=width)
  }
  return(output)
}


# Test code for convertAgevecToOutputScale
# ToDo: REMOVE THIS. THIS SHOULD NOT BE IN THE FINAL VERSION
if (FALSE){
  agevec = c(6300000,7900000,10000000,12600000,15800000,20000000,25100000,
             31600000,39800000,50100000,63100000,70800000,79400000,89100000,
             100000000,112200000,125900000,141300000,158500000,177800000,
             199500000,223900000,251200000,281800000,316200000,354800000,
             398100000,446700000,501200000,562300000,631000000,707900000,
             794300000,891300000,1000000000,1122000000,1258900000,1412500000,
             1584900000, 1778300000,1995300000,2238700000,2511900000,
             2818400000,3162300000,3548100000,3981100000,4466800000,
             5011900000,5623400000,6309600000,7079500000,7943300000,
             8912500000,10000000000,11220200000,12589300000,14125400000,
             15848900000,17782800000)
  SFR = c(4.053394, 4.056097, 4.059647, 4.064046, 4.069467, 4.076593,
          4.085262, 4.096336, 4.110347, 4.128012, 4.150411, 4.163734,
          4.178662, 4.195561, 4.214628, 4.236069, 4.260270, 4.287633,
          4.318393, 4.353160, 4.392572, 4.437295, 4.487852, 4.545175,
          4.610454, 4.684771, 4.769502, 4.866350, 4.977191, 5.104323,
          5.250973, 5.419895, 5.615865, 5.843919, 6.109948, 6.422230,
          6.790676, 7.227889, 7.750286, 8.378643, 9.140618, 10.072396,
          11.223518, 12.659710, 14.471767, 16.784642, 19.774125, 23.682655,
          28.852069, 35.743129, 44.952373, 57.121647, 72.533088, 89.625037,
          100.000000, 93.595870)
  output = convertAgevecToOutputScale(agevec, SFR, return_scale=TRUE)
  
  plot(output$age, output$data, log="x", xlim=c(4e6, 2e10), ylim=c(3, 110), pch=3, cex=2, main="log10mean")
  points(Stars$agevec, Stars$SFR, pch=4, cex=1.5)
  for (i in output$new_scale){
    abline(v=i * 1000000)
  }
}
