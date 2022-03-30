setwd("~/Documents/GitHub/LOGAN-SFH")
source("LOGAN-CLAWS.R")
library(FITSio)
library(rlang)

path_folder = "resolucion_galceran/"
waveout=seq(4700, 9400, 1.25)
colors <- c("black", "red", "blue", "green", "yellow", "pink", "cyan", "darkgreen")

FITS2Spectra <- function(filename){
  loaded <-readFITS(file.path(filename))
  wave <- seq(from=loaded$axDat$crval[1],
              to=(loaded$axDat$len[1] - 1) * loaded$axDat$cdelt[1] + loaded$axDat$crval[1],
              by=loaded$axDat$cdelt[1])
  return(list(wave=wave, flux=loaded$imDat[,1]))
}


list_files = list.files(path_folder)
list_files

originals <- c()
groupm0 <- c()
groupm1 <- c()
groupp06 <- c()
groupp15 <- c()
for (file in list_files){
  if (grepl('u.fits', file, fixed=TRUE))
    originals <- c(originals, file)
  if (grepl('m0.25', file, fixed=TRUE) & !grepl('u.fits', file, fixed=TRUE))
    groupm0 <- c(groupm0, file)
  if (grepl('m1.26', file, fixed=TRUE) & !grepl('u.fits', file, fixed=TRUE))
    groupm1 <- c(groupm1, file)
  if (grepl('p0.06', file, fixed=TRUE) & !grepl('u.fits', file, fixed=TRUE))
    groupp06 <- c(groupp06, file)
  if (grepl('p0.15', file, fixed=TRUE) & !grepl('u.fits', file, fixed=TRUE))
    groupp15 <- c(groupp15, file)
}

tmp025 = list()
tmp126 = list()
tmp006 = list()
tmp015 = list()

for (i in 1:length(groupm0)){
  tmp025[[i]] = FITS2Spectra(file.path(path_folder, groupm0[i]))
}
for (i in 1:length(groupm1)){
  tmp126[[i]] = FITS2Spectra(file.path(path_folder, groupm1[i]))
}
for (i in 1:length(groupp06)){
  tmp006[[i]] = FITS2Spectra(file.path(path_folder, groupp06[i]))
}
for (i in 1:length(groupp15)){
  tmp015[[i]] = FITS2Spectra(file.path(path_folder, groupp15[i]))
}


everything <- list(m025=list(data=tmp025,
                             original=FITS2Spectra(file.path(path_folder, "Eku1.30Zm0.25T07.5000_iTp0.00_baseFe.u.fits")),
                             name="Eku1.30Zm0.25T07.5000_iTp0.00_baseFe.u.fits",
                             list=groupm0),
                   m126=list(data=tmp126,
                             original=FITS2Spectra(file.path(path_folder, "Eku1.30Zm1.26T09.0000_iTp0.00_baseFe.u.fits")),
                             name="Eku1.30Zm1.26T09.0000_iTp0.00_baseFe.u.fits",
                             list=groupm1),
                   p006=list(data=tmp006,
                             original=FITS2Spectra(file.path(path_folder, "Eku1.30Zp0.06T14.0000_iTp0.00_baseFe.u.fits")),
                             name="Eku1.30Zp0.06T14.0000_iTp0.00_baseFe.u.fits",
                             list=groupp06),
                   p015=list(data=tmp015,
                             original=FITS2Spectra(file.path(path_folder, "Eku1.30Zp0.15T00.9000_iTp0.00_baseFe.u.fits")),
                             name="Eku1.30Zp0.15T00.9000_iTp0.00_baseFe.u.fits",
                             list=groupp15),
                   orig=originals)


data <- duplicate(everything)
difference <- duplicate(everything)

print(length(data))
for (i in 1:(length(data)-1)){
  data[[i]]$original <- interpolateToWaveout(data[[i]]$original$wave, data[[i]]$original$flux, waveout, returnList=TRUE)
  for (j in 1:length(data[[i]]$data)){
    cat(i, "-", j, "\n")
    data[[i]]$data[[j]] <- interpolateToWaveout(data[[i]]$data[[j]]$wave, data[[i]]$data[[j]]$flux, waveout, returnList=TRUE)
    difference[[i]]$data[[j]]$flux = data[[i]]$data[[j]]$flux - data[[i]]$original$flux
  }
}


# for (i in 1:(length(data)-1)) {
#   for (j in 1:length(data[[i]]$data)) {
#     cat(paste0(i, '-', j, '\n'))
#     difference[[i]]$data[[j]]$flux = data[[i]]$data[[j]]$flux - data[[i]]$original$flux
#   }
# }


for (i in 1:(length(data)-1)){
  legend_data <- c("original", data[[i]]$name)
  plot(everything[[i]]$original$wave, everything[[i]]$original$flux, type="l", main=data[[i]]$name,
       col=colors[1], ylab="Difference in flux", xlab="wavelength", xlim = c(4700, 9400))
  lines(data[[i]]$original$wave, data[[i]]$original$flux, col=colors[2])
  for (j in 1:length(data[[i]]$data)){
    lines(data[[i]]$data[[j]]$wave, data[[i]]$data[[j]]$flux, type="l", col=colors[j + 1])
    legend_data <- c(legend_data, data[[i]]$list[[j]])
  }
  legend("bottomleft", legend=legend_data, col=colors, lty=1)
  # png(file=paste0("resolucion_galceran/plots/Spectra_", data[[i]]$name, ".png"))
}


for (i in 1:(length(difference)-1)){
  yLIM <- numeric(2)
  legend_data <- c()
  for (j in 1:length(difference[[i]]$data)){
    if (max(difference[[i]]$data[[j]]$flux) > yLIM[2])
      yLIM[2] <- max(difference[[i]]$data[[j]]$flux)
    if (min(difference[[i]]$data[[j]]$flux) < yLIM[1])
      yLIM[1] <- min(difference[[i]]$data[[j]]$flux)
  }
  for (j in 1:length(difference[[i]]$data)){
    if (j == 1){
      plot(data[[i]]$data[[j]]$wave, difference[[i]]$data[[j]]$flux, main=difference[[i]]$name, col=colors[j],
           ylim=yLIM, ylab="Difference in flux", xlab="wavelength")
    } else {
      points(data[[i]]$data[[j]]$wave, difference[[i]]$data[[j]]$flux, col=colors[j])
    }
    legend_data <- c(legend_data, difference[[i]]$list[[j]])
  }
  legend("topright", legend=legend_data, col=colors, lty=1)
  png(file=paste0("resolucion_galceran/plots/Comparison_", data[[i]]$name, ".png"))
}
