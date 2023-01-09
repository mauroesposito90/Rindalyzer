## The set of custom functions created for the report module

# Load packages
library(seqinr)

## The consLen function identifies the length of the consensus sequence
# @param consRef Link to consensus reference
# @return length Length of the consensus sequence
#
consLen <- function(consRef){
    fasta <- read.fasta(file = consRef)
    length <- getLength(object = fasta)
    return(length)
}

## The logUpdater function create/update the log file
# @param modality Create or Update the log file
# @param message Message to write in the log file
# @param outDir Directory for the output
#
logUpdater <- function(modality, message, outDir){
    if(modality=="Create"){
        file.create("Rindalyzer.log")
        appendOption=FALSE
    }else if(modality=="Update"){
        appendOption=TRUE
    }
    write(x = message, file = paste0(outDir, "/Rindalyzer.log"), 
          append = appendOption)
}