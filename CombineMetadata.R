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
for (arg in args){
    if (substr(arg, 1, 4) == "out="){
        outfile = substr(arg, 5, nchar(arg))
    } else if (substr(arg, 1, 5) == "path="){
        filepath = substr(arg, 6, nchar(arg))
    }
}


if (is.null(outfile)){
    cat("No value for output filename (out=) was given. Output filename set to default (MetaD_combined.rda).\n")
} else{
    cat(paste0("Output filename = ", outfile, "\n"))
    outfile="MetaD_combined.rda"
}

# ToDo: Fix this hardcoded stuff... environ variables??
if (is.null(filepath)){
    cat("No value for filepath (path=) was given. Filepath set to default (/Users/enrique/Documents/GitHub/LOGAN-SFH/).\n")
} else{
    cat(paste0("Output filename = ", "/Users/enrique/Documents/GitHub/LOGAN-SFH/", "\n"))
    filepath="/Users/enrique/Documents/GitHub/LOGAN-SFH/"
}



metadata_to_save = list()








