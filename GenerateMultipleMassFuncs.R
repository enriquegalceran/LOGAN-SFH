##############
# Generate multiple spectra based on a specific massfunction.

library(ProSpect)
library(stringi)
library(jsonlite)
source("GenerateOutputsFinal.R")


generateSpecFromParams <- function(massParams="default",
                                   folderPath="OutputsGenerated",
                                   ZParams="default",
                                   speclib=NULL,
                                   stellpop="EMILESCombined",
                                   filters="default",
                                   emission=TRUE,
                                   emission_scale="SFR",
                                   waveout=seq(4700, 9400, 1.25),
                                   CRVALDELTA="default",
                                   forcemass=1e10,
                                   notesParams=NULL,
                                   randomSamples=0,
                                   SNRatio=30,
                                   onlyNoise=FALSE,
                                   addDate=TRUE,
                                   filenamemode=2,
                                   confirmation=FALSE,
                                   verbose=2,
                                   verboseSteps=1,        # if 0, will be the same as pythonStep
                                   cleanOutputFolder=FALSE,
                                   absolutePath=FALSE,
                                   bytesForPython=50e6,
                                   singleOutput=TRUE,
                                   metadataName="MetadataOutput.json",
                                   ...) {
  #####
  # Functions
  advanceOneValue <- function(argIdx){
    n = length(argIdx$idx)
    argIdx$idx[n] = argIdx$idx[n] + 1
    i = n
    while (argIdx$idx[i] > argIdx$maxIdx[i] && i > 1){
      argIdx$idx[i] = 1
      i = i - 1
      argIdx$idx[i] = argIdx$idx[i] + 1
    }
    return(argIdx)
  }
  
  
  interpolateToWaveout <- function(inputMatrix, waveout){
    waveoutL = log10(waveout)
    x1L = log10(inputMatrix[,1])
    y1L = log10(inputMatrix[,2])
    spect = 10^approxfun(x1L, y1L, rule=2)(waveoutL)
    return(data.frame(wave=waveout, flux=spect))
  }
  
  
  argumentsToHeaders <- function(params, noteParams=NULL){
    out <- list(keyword=c(), value=c(), note=rep("", length(params)))
    for (i in 1:length(params)){
      out$keyword = c(out$keyword, names(params)[i])
      out$value = c(out$value, params[[i]])
      if (!is.null(noteParams)){
        if (names(params)[i] %in% names(noteParams)){
          out$note[i] = noteParams[[names(params)[i]]]
        }
      }
    }
    return(out)
  }
  
  
  generateFilename <- function(params, mode=2, rngLength=8, addDate=TRUE, singleOutput=FALSE){
    if (mode == 1){
      filename = paste0(params$mfunc, "_")
      if (addDate){
        now = Sys.time()
        filename = paste0(filename, gsub(" ", "T", now), "_")
      }
      
      ZArgs = c("Zstart", "Zfinal", "yield", "Zagemax", "mfunc", "rndSample")
      nonZArgs = names(params)[!(names(params) %in% ZArgs)]
      # Add mfunc params
      for (i in 1:length(nonZArgs)){
        if (params[[nonZArgs[i]]]>0){
          union = "p"
        } else {
          union = "m"
        }
        filename = paste(filename, abs(params[[nonZArgs[i]]]), sep=union)
      }
      # Add Z params
      filename = paste0(filename, "_")
      for (i in 1:4){
        if (params[[ZArgs[i]]]>0){
          union = "p"
        } else {
          union = "m"
        }
        filename = paste(filename, abs(params[[ZArgs[i]]]), sep=union)
      }
      # Add idx for random sample
      filename = paste0(filename, "_", params[["rndSample"]])
      
      # Verify no + or - symbols are left
      filename = gsub("\\+", "p", filename)
      filename = gsub("\\-", "m", filename)
      
    } else if (mode == 2){
      if (addDate){
        now = Sys.time()
        if (!singleOutput){
          filename = paste(params$mfunc, gsub(" ", "T", gsub(":", "", gsub("-", "", now))), stri_rand_strings(n=1, length=rngLength), sep="_")
        } else {
          filename = paste(gsub(" ", "T", gsub(":", "", gsub("-", "", now))), stri_rand_strings(n=1, length=rngLength), sep="_")
        }
      } else {
        filename = paste(params$mfunc, stri_rand_strings(n=1, length=rngLength), sep="_")
      }
    } else if (mode == 3){
      stop("mode==3 -> Aun no implementado!!!")
    } else if (mode == 4){
      filename = stri_rand_strings(n=1, length=16)
    }
    return(filename)
  }
  
  insertNoise <- function(spec, SNRatio){
    # Each point in the spectrum will be modified
    # according to it's value and the SNRatio.
    spec$flux$flux = rnorm(length(spec$flux$flux), spec$flux$flux, spec$flux$flux/SNRatio)
    return(spec)
  }
  
  roundMultiple <- function(number, multiple, f=round){
    return(f(number/multiple)*multiple)
  }
  
  bytes2Human <- function(B, decimals=2){
    betterUnits = roundMultiple(log10(B), 3, trunc)
    B = B/(10^betterUnits)
    units = c("", "K", "M", "G", "T", "P")
    units = units[betterUnits/3 + 1]
    return(paste0(toString(round(B, decimals)), units, "B"))
  }
   
   
  #####
  # Preprocess data
  dots = list(...)
  
  # Select spectra library
  if (!is.null(speclib)){
    if (!is.null(stellpop) & stellpop != "EMILESCombined"){
      warning("stellpop and speclib do not match!")
    } else {
      # Return to one of the default accepted values for SFHfunc
      stellpop="EMILES"
    }
  } else {
    if (stellpop  == "EMILESCombined"){
      # Return to one of the default accepted values for SFHfunc
      stellpop = "EMILES"
      # TODO: This will need to read global variable/env-variable or something similar.
      speclib = readRDS(file="EMILESCombined.rds")
    }
  }
  
  # Select filter data
  if (filters == "default"){
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
  }
  
  # Select massParams
  if (is.character(massParams)){
    if (massParams == "default"){
      massParams = list(
        dtau=list(
          name="dtau",
          func=massfunc_dtau
        )
      )
    } else {
      print(massParams)
      stop("ERROR IN massParams!^")
    }
  }
  
  # Select ZParams
  if (is.character(ZParams)){
    if (ZParams == "default"){
      ZParams = list(
        func=Zfunc_massmap_box,
        Zstart=1e-4,
        Zfinal=0.02,
        yield=0.03
      )
    } else {
      print(ZParams)
      stop("ERROR IN ZParams!^")
    }
  } else {
    if (any(sort(names(ZParams)) != sort(c("func", "Zstart", "Zfinal", "yield")))){
      stop("ZParams' arguments need to be: func, Zstart, Zfinal and yield")
    }
    if (any(names(ZParams) != c("func", "Zstart", "Zfinal", "yield"))){
      # Rearrange parameters
      ZParams = ZParams[c("func", "Zstart", "Zfinal", "yield")]
    }
  }
  
  # AbsolutePath
  if (absolutePath){
    folderPath = folderPath
  } else {
    folderPath = file.path(getwd(), folderPath)
  }
  if (verbose > 0)
    cat(paste0("Output Folder: ", folderPath, "\n"))
  
  # Verify if folder exists (and create it if it does not)
  if (!dir.exists(folderPath)){
    dir.create(folderPath)
  }
  
  nfilesAtStart <- length(list.files(folderPath))
  
  # Clean Output Folder
  if (cleanOutputFolder & nfilesAtStart > 0){
    if (confirmation){
      cleanConfirmation <- readline(prompt = paste0("Clean Outputfolder with ", nfilesAtStart, " elements? (Y/[N]/All) "))
      if (any(cleanConfirmation == c("y", "Y", "yes", "YES", "all", "All", "ALL", "a"))){
        cleanConfirmation = TRUE
      } else {
        cleanConfirmation = FALSE
      }
    } else {
      cleanConfirmation = TRUE
    }

    if (cleanConfirmation){
      ls = list.files(folderPath)
      for (file in ls){
        if (any(cleanConfirmation == c("all", "ALL", "All", "a")) || (file != "ReferenceInputs.fits" & file != "ReferenceLabel.fits"))
          file.remove(file.path(folderPath, file))
      }
      if (verbose > 0)
        cat("Output Folder has been cleaned.\n")
    }
  }
  
  # Verify if the Reference files for Inputs and Labels is available.
  ls = list.files(folderPath)
  if (!("ReferenceInput.fits" %in% ls)){
    # Generate Reference Input
    wave = data.matrix(seq(4700, 9400, 1.25))
    writeFITSim(wave,
                file = file.path(folderPath, "ReferenceInput.fits"),
                ctypen = c("wave", ""),
                cunitn = c("A", ""),
                crvaln = c(4700, 1),
                cdeltn = c(1.25, 1),
                c1 = "Reference Inputs. Contains wavelength (see crval&cdelt)")
    cat("Reference file for Inputs generated.\n")
  }
    
  agevec = speclib$Age
  if (!("ReferenceLabel.fits" %in% ls)){
    # Generate Reference Input
    # agevec = data.matrix(c(6300000,7900000,10000000,12600000,15800000,20000000,25100000,
    #                        31600000,39800000,50100000,63100000,70800000,79400000,89100000,
    #                        100000000,112200000,125900000,141300000,158500000,177800000,
    #                        199500000,223900000,251200000,281800000,316200000,354800000,
    #                        398100000,446700000,501200000,562300000,631000000,707900000,
    #                        794300000,891300000,1000000000,1122000000,1258900000,1412500000,
    #                        1584900000,1778300000,1995300000,2238700000,2511900000,
    #                        2818400000,3162300000,3548100000,3981100000,4466800000,
    #                        5011900000,5623400000,6309600000,7079500000,7943300000,
    #                        8912500000,10000000000,11220200000,12589300000,14125400000,
    #                        15848900000,17782800000))
    writeFITSim(agevec,
                file = file.path(folderPath, "ReferenceLabel.fits"),
                ctypen = c("Agevec", ""),
                cunitn = c("yrs", ""),
                c1 = "Reference Labels. Contains Ages (non-linear)")
    cat("Reference file for Labels generated.\n")
  }
  
  # Print which parameters are going to be edited
  if (verbose > 0)
    cat("\n -------------- Parameters that will change  -------------- \n")
  totalNumberOfCases = 0
  
  if (verbose > 0)
    cat("Parameters for Z:\n")
  combZ = 1
  for (zPar in 2:4){
    if (verbose > 0){
      if (length(ZParams[[zPar]]) > 1){
        cat(sprintf("%s: [%s,%s] (%d)\n",
                    names(ZParams)[zPar],
                    ZParams[[zPar]][1],
                    ZParams[[zPar]][length(ZParams[[zPar]])],
                    length(ZParams[[zPar]])))
      } else {
        cat(sprintf("%s: [%s] (%d)\n",
                    names(ZParams)[zPar],
                    ZParams[[zPar]][1],
                    length(ZParams[[zPar]])))
      }
    }
    combZ = combZ * length(ZParams[[zPar]])
  }
  if (verbose > 0)
    cat(sprintf("Total combinations for Z: %d\n ----------- \n", combZ))
  
  for (funct in massParams){
    if (verbose > 0)
      cat(sprintf("Function: %s\n", funct$name))
    combMassFunc = 1
    if (length(funct)>2){
      for (param in 3:length(funct)){
        if (verbose > 0){
          if (length(funct[[param]])>1){
            cat(sprintf("%s: [%s,%s] (%d)\n",
                        names(funct)[param],
                        funct[[param]][1],
                        funct[[param]][length(funct[[param]])],
                        length(funct[[param]])))
          } else {
            cat(sprintf("%s: [%s] (%d)\n",
                        names(funct)[param],
                        funct[[param]][1],
                        length(funct[[param]])))
          }
        }
        combMassFunc = combMassFunc * length(funct[[param]])
      }
    }
    if (verbose > 0)
      cat(sprintf("Total combinations for %s: %d\n", funct$name, combMassFunc))
    totalNumberOfCases = totalNumberOfCases + combMassFunc
    if (verbose > 0)
      cat(" ----------- \n")
  }
  totalNumberOfCases = totalNumberOfCases*combZ*(randomSamples + 1)
  if (verbose > 0) {
    cat(paste0("Number of Resamplings with noise: ", randomSamples, " (+1[no noise] = ", randomSamples + 1, ")\n"))
    cat(sprintf("Total number of cases (massfuncs + Zs + rndSamples): %d\n", totalNumberOfCases))
    if (CRVALDELTA == "default"){
      expectedSizePerPair = 43200
    } else {
      expectedSizePerPair = 72000
    }
    cat(sprintf(paste0("Approximate expected total file size (default filters&waveout): ", bytes2Human(totalNumberOfCases*expectedSizePerPair), "\n\n")))
  }
 
   
  #####
  # Calculate how many pairs are to be generated before calling the Python code.
  pythonSteps = trunc(bytesForPython/expectedSizePerPair)
  cat("The Python script will be executed every ", pythonSteps, " calculations. (", bytes2Human(expectedSizePerPair * pythonSteps), "/", bytes2Human(bytesForPython), ")\n", sep="")
  pythonListFilenames = rep(NA, min(pythonSteps, totalNumberOfCases))
  if (verboseSteps == 0){
    verboseSteps = pythonSteps
  }
  
  
  if (confirmation){
    confirmation <- readline(prompt = paste0("Continue operation for ", totalNumberOfCases, " combinations? ([Y]/N) "))
    if (!any(confirmation == c("", "y", "Y", "yes", "YES"))){
      stop("User stopped execution.")
    }
  }
  
  # Initiate the matrix that will contain all of the data
  numberColumnsIn = length(waveout) + length(filters) + 1    # Add 1 for the ID
  # numberColumnsLa = 2 * length(agevec) + 1                   # Add 1 for the ID
  numberColumnsLa = 2 * 56 + 1                               # 56 should be the length of the Agevector, but it includes very long sizes, so this may generate some problems in the future TODO check this!
  completeDataMatrixIn <- matrix(, nrow=totalNumberOfCases + 1, ncol= numberColumnsIn)
  completeDataMatrixLa <- matrix(, nrow=totalNumberOfCases + 1, ncol= numberColumnsLa)
  # ToDo: set the first row correctly (Corregir lo de 56)
  completeDataMatrixIn[1, ] = c(0, waveout, seq(1:length(filters)))
  completeDataMatrixLa[1, ] = c(0, agevec[1:56], agevec[1:56])
  
  # Order of parameters for the metadata
  orderParameters = c()
  filterData = NULL
  
  #####
  # Iterate over massParams
  absoluteCountCases = 1
  # Time the code from here (after cleaning/not cleaning the folder)
  ptm <- proc.time()
  for (func in massParams){
    if (verbose > 0)
      cat(paste0("\n -------------- Starting to calculate for massfunction: ", func$name, " -------------- \n"))
    # Verify if there are additional parameters
    if (length(func) > 2){
      massMoreThan1 = lengths(func[3:length(func)])>1
      # Separate into single element and multiple elements
      massMultElements = func[names(massMoreThan1)[massMoreThan1]]
      massOneElement = func[names(massMoreThan1)[!massMoreThan1]]
    } else {
      massMultElements = list()
      massOneElement = list()
    }
    
    # Generate a list of all default parameters that have not been set
    massfuncDefaultArgs <- formals(func$func)[!(names(formals(func$func)) %in% names(func))]
    # Remove Symbols
    massfuncDefaultArgs <- massfuncDefaultArgs[!sapply(massfuncDefaultArgs, is.symbol)]
    
    # For parameters that have multiple values, we will iterate over them.
    massArgIdx = list(paramName=names(massMultElements), idx = rep(1, length(massMultElements)), maxIdx = lengths(massMultElements))
    massNumberOfCases = prod(massArgIdx$maxIdx)
    if (massNumberOfCases == 0)
      massNumberOfCases = 1
    massCaseNumber = 1L
    
    
    #####
    # Iterate over all the possible elements
    while (massCaseNumber <= massNumberOfCases){
      
      # Generate the full list of arguments
      # Combine default elements + static elements
      massfunc_args = c(massfuncDefaultArgs, massOneElement)
      # Evaluate current indexes
      massCurrentArgs = list()
      if (length(massArgIdx$idx)>0) {
        for (k in 1:length(massArgIdx$idx)){
          massCurrentArgs = c(massCurrentArgs, list(massMultElements[[k]][massArgIdx$idx[[k]]]))
        }
      }
      names(massCurrentArgs) <- massArgIdx$paramName
      # Combine with previous values
      massfunc_args <- c(massfunc_args, massCurrentArgs)
      # Sort according to default values
      massfunc_args <- massfunc_args[names(formals(func$func))[names(formals(func$func)) %in% names(massfunc_args)]]
      
      
      #####
      # Calculating Arguments for Z
      if (length(ZParams) > 1){
        ZMoreThan1 = lengths(ZParams[2:length(ZParams)])>1
        # Separate into single element and multiple elements
        ZMultElements = ZParams[names(ZMoreThan1)[ZMoreThan1]]
        ZOneElement = ZParams[names(ZMoreThan1)[!ZMoreThan1]]
      } else {
        ZMultElements = list()
        ZOneElement = list()
      }
      
      # Generate a list of all default parameters that have not been set
      ZfuncDefaultArgs <- formals(ZParams$func)[!(names(formals(ZParams$func)) %in% names(ZParams))]
      # Remove Symbols
      ZfuncDefaultArgs <- ZfuncDefaultArgs[!sapply(ZfuncDefaultArgs, is.symbol)]
      
      # For parameters that have multiple values, we will iterate over them.
      ZArgIdx = list(paramName=names(ZMultElements), idx = rep(1, length(ZMultElements)), maxIdx = lengths(ZMultElements))
      ZNumberOfCases = prod(ZArgIdx$maxIdx)
      if (ZNumberOfCases == 0)
        ZNumberOfCases = 1
      ZCaseNumber = 1L
      
      
      #####
      # Iterate over all the ZParameters defined
      while (ZCaseNumber <= ZNumberOfCases){
        # Generate the full list of arguments
        # Combine default elements + static elements
        Zfunc_args = c(ZfuncDefaultArgs, ZOneElement)
        # Evaluate current indexes
        ZCurrentArgs = list()
        if (length(ZArgIdx$idx)>0) {
          for (k in 1:length(ZArgIdx$idx)){
            ZCurrentArgs = c(ZCurrentArgs, list(ZMultElements[[k]][ZArgIdx$idx[[k]]]))
          }
        }
        names(ZCurrentArgs) <- ZArgIdx$paramName
        # Combine with previous values
        Zfunc_args <- c(Zfunc_args, ZCurrentArgs)
        # Sort according to default values
        Zfunc_args <- Zfunc_args[names(formals(ZParams$func))[names(formals(ZParams$func)) %in% names(Zfunc_args)]]


        
        #####
        # Evaluate and calculate SFH function
        spectraObject = do.call('SFHfunc', c(list(massfunc = func$func,
                                                  forcemass = forcemass,
                                                  stellpop = stellpop,
                                                  speclib = speclib,
                                                  filters = filters,
                                                  emission = emission,
                                                  emission_scale = emission_scale
                                                  ),
                                             massfunc_args,
                                             Zfunc_args
                                             )
                                )
        
        
        #####
        # Generate new data according to errors defined
        # Error resampling is done BEFORE waveout
        for (rnd in 0:randomSamples){
          if (rnd > 0){
            spectraObject = insertNoise(spectraObject, SNRatio)
          }
        
          
          #####
          # Adjust spectra's Wavelength according to waveout if not NULL
          if (!is.null(waveout)){
            spectraObject$flux = interpolateToWaveout(spectraObject$flux, waveout)
          } 
          

          #####
          # Export file
          # First calculate the filename that is going to be used
          savedParams = c(list(mfunc=func$name), massfunc_args, Zfunc_args, list(rndSample=rnd))
          filename = generateFilename(savedParams, mode=filenamemode, rngLength=6, addDate = addDate)
          notesParams = list(mSFR="test1234", yield="Another test note with spaces", theforcebewithyou="keyword that is not inthe other params")
          savedParams = argumentsToHeaders(savedParams, notesParams)
          
          # Output through screen where the code is at.
          # Calculate length of strings
          ncharsNumberCases = sapply(sapply(c(massNumberOfCases, ZNumberOfCases, randomSamples + 1), toString), nchar)
          sFilename = sprintf("%s - ", filename)
          sMass = sprintf(paste0("mF:(%0", ncharsNumberCases[1],"d/%d) - "), massCaseNumber, massNumberOfCases)
          sZ = sprintf(paste0("Z:(%0", ncharsNumberCases[2],"d/%d) - "), ZCaseNumber, ZNumberOfCases)
          sRNG = sprintf(paste0("RNG:(%0", ncharsNumberCases[3],"d/%d) - "), rnd + 1, randomSamples + 1)
          sCurrentMF = sprintf(paste0("C.mF:(%0", nchar(toString(massNumberOfCases * ZNumberOfCases * (randomSamples + 1))), "d/%d) - "),
                           (rnd + 1) + (ZCaseNumber - 1) * (randomSamples + 1) + (massCaseNumber - 1) * ZNumberOfCases * (randomSamples + 1),
                          massNumberOfCases * ZNumberOfCases * (randomSamples + 1)
                          )
          sTotal = sprintf(paste0("Tot:(%0", nchar(toString(combMassFunc*combZ*(randomSamples + 1))), "d/%d) - %5.1f"),
                           absoluteCountCases,
                           totalNumberOfCases,
                           absoluteCountCases/totalNumberOfCases*100
                           )
          if (verbose > 0){
            if (absoluteCountCases %% verboseSteps == 0 || absoluteCountCases == 1 || absoluteCountCases == totalNumberOfCases){
              cat(paste0(sFilename, sMass, sZ, sRNG, sCurrentMF, sTotal, "%\n"))
            }
          }
          
          # EXPORT FILE
          if (!singleOutput){
            # Generate a file for each element
            exportObjectToFITS(spectraObject,
                               filename = filename,
                               foldername = folderPath,
                               spectrumParam = savedParams,
                               randomNoise = c(rnd, SNRatio),
                               verbose = verbose,
                               absolutePath = TRUE,
                               CRVALDELTA = CRVALDELTA,
                               forcemass = forcemass
                               )
          } else {
            # Generate a single for for every element
            # ToDo: Generate a function that, given an ID and parameters, gives out the parameters and generates a single file based on the data provided
            
            # New row to be added (Input)
            # ID, spectra, Magnitudes
            newRowIn <- c(absoluteCountCases, spectraObject$flux$flux, spectraObject$out$out)  # ESTO ES LO QUE SALE NULL!!! # TODO
            newRowLa <- c(absoluteCountCases, spectraObject$SFR, spectraObject$Zvec)
            
            # Add new row to the Matrix
            completeDataMatrixIn[absoluteCountCases + 1, ] = newRowIn
            completeDataMatrixLa[absoluteCountCases + 1, ] = newRowLa
          }
  
          # ToDo: Fix every part of the code that depends on "singleOutput"
          pythonListFilenames[absoluteCountCases %% pythonSteps] = filename
          #####
          # Execute Python code that reduces to 32-bit every pythonStep
          if ((absoluteCountCases %% pythonSteps == 0 | absoluteCountCases == totalNumberOfCases) && !singleOutput){
            listFilename = file.path(getwd(), "ListToConvertTo32.txt")
            zz <- file(listFilename, "wb")
            writeBin(paste(pythonListFilenames, collapse="\n"), zz)
            close(zz)
            # TODO: output something to screen IF VERBOSE
            pythonExecution = paste0("/Users/enrique/.conda/envs/LOGAN/bin/python ConsolidateTables.py --list ", listFilename, " --datadirectory ", folderPath)
            if (verbose>0)
              cat("Executing: '", pythonExecution, "'\n", sep="")
            system(pythonExecution)
            pythonListFilenames = c()
            # TODO: Using env variables There should be an env variable pointing to the python execution
            # TODO: Should be more flexible
          }
          
          # Count every case
          absoluteCountCases = absoluteCountCases + 1
        }
        ZArgIdx <- advanceOneValue(ZArgIdx)
        ZCaseNumber <- ZCaseNumber + 1L
      }
      massArgIdx <- advanceOneValue(massArgIdx)
      massCaseNumber <- massCaseNumber + 1L
    }
    if (verbose > 0)
      cat(paste0(" -------------- Finished to calculate for massfunction: ", func$name, " -------------- \n"))
    
    # Add information to orderParameters (for metadata file)
    tmp = list(list(mass = massArgIdx$paramName, Z = ZArgIdx$paramName))
    names(tmp) <- func$name
    orderParameters = c(orderParameters, tmp)
  }
  
  #########
  # Store final Matrix into a FITS "image" in 2D
  if (singleOutput){
    filename = generateFilename(NULL, mode=filenamemode, rngLength=6, addDate = addDate, singleOutput=TRUE)
    
    filterData = spectraObject$out
    UUIDs <-exportObjectsToSingleFITS(inputMatrix = completeDataMatrixIn,
                                      labelMatrix = completeDataMatrixLa,
                                      filename = filename,
                                      foldername = folderPath,
                                      filters = filterData,
                                      absolutePath = absolutePath,
                                      verbose=verbose
    )
    # ToDo: Reduce data with Python
    # ToDo: Make Python understand that there are two distinct functions and act accordingly
    
    # Export configuration settings into a JSON file
    jsonData <- toJSON(list(massParams=massParams,
                            ZParams=ZParams,
                            orderParameters = orderParameters,
                            filters=names(filters),
                            emission=emission,
                            emission_scale=emission_scale,
                            forcemass=forcemass,
                            randomSamples=randomSamples,
                            SNRatio=SNRatio,
                            onlyNoise=onlyNoise,
                            UUIDs = UUIDs
                            )
                       )
    write(jsonData, metadataName)
    if (verbose >= 1)
      cat("Metadata was stored in", metadataName, "...\n")
    
  }
  
  if (verbose > 0){
    cat("\n\n")
    filesAtEnd = length(list.files(folderPath))
    cat(paste0("New Files Generated: ",  filesAtEnd - nfilesAtStart, "\n"))
    cat(paste0("Individual Cases Computed: ",  (filesAtEnd - nfilesAtStart) / 2, "\n"))
    cat("Time:\n")
    print(proc.time() - ptm)
    cat(" -------- FINISHED --------\n")
  }
  
}
