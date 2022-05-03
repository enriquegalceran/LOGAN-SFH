#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==1) {
    # default output file
    args[2] = "out=MetaD_combined.rda"
}

# Select outfile and filepath
outfile=NULL
filepath=NULL
uuidi=NULL
uuidl=NULL
uuidm=NULL
for (arg in args){
    if (substr(arg, 1, 4) == "out="){
        outfile = substr(arg, 5, nchar(arg))
    } else if (substr(arg, 1, 5) == "path="){
        filepath = substr(arg, 6, nchar(arg))
    } else if (substr(arg, 1, 6) == "uuidi="){
        uuidi = substr(arg, 7, nchar(arg))
    } else if (substr(arg, 1, 6) == "uuidl="){
        uuidl = substr(arg, 7, nchar(arg))
    } else if (substr(arg, 1, 6) == "uuidm="){
        uuidm = substr(arg, 7, nchar(arg))
    } 
}
for (arg in args){
    if (substr(arg, 1, 4) == "uuid")
        args = args[-which(args==arg)]
}


if (is.null(outfile)){
    cat("No value for output filename (out=) was given. Output filename set to default (MetaD_combined.rda).\n")
    outfile="MetaD_combined.rda"
} else{
    args = args[-which(args == paste0("out=", outfile))]
}
cat(paste0("Output file = ", outfile, "\n"))


# TODO: Fix this hardcoded stuff... environ variables??
if (is.null(filepath)){
    cat("No value for filepath (path=) was given. Filepath set to current directory (", getwd(), ").\n")
    filepath=getwd()
} else{
    args = args[-which(args == paste0("path=", filepath))]
    if (substr(filepath, 1, 1) != "/")
        filepath = file.path(getwd(), filepath)
    cat(paste0("Output filename = ", filepath, "\n"))
}

# Add path to outfile
if (substr(outfile, 1, 1) != "/"){
    outfile = file.path(filepath, outfile)
}

metadata_to_save = list()
metadata_to_save = append(metadata_to_save, list("-1"=list(uuidi=uuidi, uuidl=uuidl, uuidm=uuidm)))

i = 0
for (arg in args){
    
    # inputname = file.path(filepath, paste0("Input_", args, ".fits"))
    # labelname = file.path(filepath, paste0("Label_", args, ".fits"))
    MetaDname = file.path(filepath, paste0("MetaD_", arg, ".rda"))
    
    cat(paste0("Loading file '", MetaDname, "' ...\n"))
    load(MetaDname)
    
    metadata_to_save = append(metadata_to_save, list(metadata))
    
    i = i + 1
    
}
metadata_combined = setNames(metadata_to_save, seq(0, (i-1)))
cat(paste0("Saving output in ", outfile, " ...\n"))
save(metadata_combined, file=outfile)

