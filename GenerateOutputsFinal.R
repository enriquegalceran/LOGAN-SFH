##############
# This code generates from an object the expected output files
# These files are 
#       1) a .fits (INPUT) file with a table which includes the spectra and the fotometry in the header
#       2) a .fits (LABEL) file with the SFR generated (Solution)
# First, the data is recollected, then the additional headers are created, the main header is generated, and finally the file is saved.

# extraHeader format: list with elements: "keyword", "value", "note"
# spectrumParam format: idem extraHeader, exclusive for parameters for Spectrum

#####
# Load required libraries
library(FITSio)
library(ProSpect)
library(uuid)

exportObjectToFITS <- function(inputObject,
                               filename,
                               foldername,
                               spectrumParam,
                               extraHeaders = NULL,
                               randomNoise = c(0, 0),
                               fileprefix = "",
                               absolutePath = FALSE,
                               verbose=2,
                               CRVALDELTA="default",
                               forcemass=FALSE
                               ) {
  #####
  # Function(s)
  # Adds the Filters to the header
  addFiltersToHeader <- function(header, filters){
    header <- addComment("------ List of filters analysed ------", header=header)
    for (i in 1:length(filters$filter)){
      header <- addKwv(paste0("filterN", i),
                       filters$filter[i],
                       note=paste("Filter Name", i),
                       header=header)
      header <- addKwv(paste0("filterV", i),
                       filters$out[i],
                       note=paste("Filter Value", i),
                       header=header)
      header <- addKwv(paste0("filterC", i),
                       filters$cenwave[i],
                       note=paste("Filter CenterWave", i),
                       header=header)
    }
    return(header)
  }
  
  addExtraHeaderElements <- function(spectrumParam=NULL, extraHeaders, header){
    if (!is.null(spectrumParam)){
      header <- addComment("------ MassFunction details ------", header=header)
      for (i in 1:length(spectrumParam$keyword)){
        if (suppressWarnings(!is.na(as.numeric(spectrumParam$value[i])))){
          header <- addKwv(spectrumParam$keyword[i],
                           suppressWarnings(as.numeric(spectrumParam$value[i])),
                           spectrumParam$note[i],
                           header=header)
        } else {
          header <- addKwv(spectrumParam$keyword[i],
                           spectrumParam$value[i],
                           spectrumParam$note[i],
                           header=header)
        }
      }
    }
    if (!is.null(extraHeaders)){
      for (i in 1:length(extraHeaders)){
        if (suppressWarnings(!is.na(as.numeric(extraHeaders$value[i])))){
          header <- addKwv(extraHeaders$keyword[i],
                           suppressWarnings(as.numeric(extraHeaders$value[i])),
                           extraHeaders$note[i],
                           header=header)
        } else {
          header <- addKwv(extraHeaders$keyword[i],
                           extraHeaders$value[i],
                           extraHeaders$note[i],
                           header=header)
        }
      }
    }
    return(header)
  }

  #####
  # Set default values
  if (CRVALDELTA == "default"){
    crvalI1 = 4700
    cdeltI1 = 1.25
  } else {
    crvalI1 = 1
    cdeltI1 = 1
  }
  
  
  #####
  # Read and generate data + metadata
  if (CRVALDELTA == "default"){
    spectra = data.matrix(inputObject$flux$flux)
    Labels = cbind(inputObject$SFR, inputObject$Zvec)
  } else {
    spectra = data.matrix(inputObject$flux)
    Labels = cbind(inputObject$agevec, inputObject$SFR, inputObject$Zvec)
  }
  filters = inputObject$out
  objectUUID = UUIDgenerate(n=2L)
  dateGenerated = Sys.time()
  
  
  #####
  # Generate Header for Inputs
  hdrIn <- newKwv("DateGen", dateGenerated, "Date when generated")
  hdrIn <- addKwv("FileType", "Input", note="File type (Input/Label)", header=hdrIn)
  hdrIn <- addKwv("Filename", paste0("Input_", fileprefix, filename, ".fits"), header=hdrIn)
  hdrIn <- addKwv("UuidInp", objectUUID[1], note="UUID for Input", header=hdrIn)
  hdrIn <- addKwv("UuidLab", objectUUID[2], note="UUID for Label", header=hdrIn)
  hdrIn <- addExtraHeaderElements(spectrumParam, extraHeaders, hdrIn)
  hdrIn <- addKwv("Minwave", min(spectra[,1]), note="Minimum wavelength", header=hdrIn)
  hdrIn <- addKwv("Maxwave", max(spectra[,1]), note="Maximum wavelength", header=hdrIn)
  hdrIn <- addKwv("Nelement", length(spectra[,1]), note="Elements in wavelength", header=hdrIn)
  hdrIn <- addKwv("masstot", inputObject$masstot, note="Total mass", header=hdrIn)
  hdrIn <- addKwv("forcema", if(forcemass==FALSE){FALSE} else {TRUE}, note="Was the mass forced to a value", header=hdrIn)
  hdrIn <- addKwv("Noise", if(randomNoise[1] >0){1}else{0}, note="0: No noise  - 1: Noise according to SNR", header=hdrIn)
  hdrIn <- addKwv("SNR", randomNoise[2], note="Signal-Noise Ratio", header=hdrIn)
  hdrIn <- addFiltersToHeader(hdrIn, filters)
  
  
  #####
  # Generate Header for Labels
  hdrLb <- newKwv("DateGen", dateGenerated, "Date when generated")
  hdrLb <- addKwv("FileType", "Label", note="File type [Input/Label]", header=hdrLb)
  hdrLb <- addKwv("Filename", paste0("Label_", fileprefix, filename, ".fits"), header=hdrLb)
  hdrLb <- addKwv("UuidInp", objectUUID[1], note="UUID for Input", header=hdrLb)
  hdrLb <- addKwv("UuidLab", objectUUID[2], note="UUID for Label", header=hdrLb)
  hdrLb <- addComment("--- Column Values ---", header=hdrLb)
  hdrLb <- addKwv("TTYPE1", "Agevector", header=hdrLb)
  hdrLb <- addKwv("TTYPE2", "SFR", header=hdrLb)
  hdrLb <- addKwv("TTYPE3", "Z", header=hdrLb)
  hdrLb <- addExtraHeaderElements(extraHeaders=extraHeaders, header=hdrLb)
  
  
  #####
  # Generate File Names and verify directory
  if (absolutePath){
    filedirectory = foldername
  } else {
    filedirectory = file.path(getwd(), foldername)
  }
  filename = paste0(fileprefix, filename, ".fits")
  
  # If output folder does not exist, generate it
  dir.create(filedirectory, showWarnings = FALSE)
  if (verbose > 1)
    cat("Saving files in '", filedirectory, "/'.\n", sep="")
  
  # New file name system does not use this. If error appears with filename,
  #     this should be uncommented, having added previously a hidewarnings/similar
  # # Check if files exist and remove if necessary
  # if (file.exists(paste0(filedirectory, "/Input_", filename))) {
  #   #Delete file if it exists
  #   file.remove(paste0(filedirectory, "/Input_", filename))
  # }
  # if (file.exists(paste0(filedirectory, "/Label_", filename))) {
  #   #Delete file if it exists
  #   file.remove(paste0(filedirectory, "/Label_", filename))
  # }
  
  
  #####
  # Save files
  if (verbose > 1)
    cat(paste0("Saving Input_", filename, " ...\n"))
  writeFITSim(spectra,
              file = paste0(filedirectory, "/Input_", filename),
              ctypen = c("wave", "flux"),
              cunitn = c("A", "erg/s/cm**2/A"),
              crvaln = c(crvalI1, 1),
              cdeltn = c(cdeltI1, 1),
              header=hdrIn)
  
  if (verbose > 1)
    cat(paste0("Saving Label_", filename, " ...\n"))
  writeFITSim(Labels,
              file = paste0(filedirectory, "/Label_", filename),
              ctypen = c("agevec", "SFR/Z"),
              cunitn = c("Gyrs", "Msol/yrs"),
              header=hdrLb)
}
