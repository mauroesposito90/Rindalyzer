## The function of the module fastq is to create the FASTQ from the consensus BAM
# @param consensusBAM Path to the consensus BAM
# @param nThreads Number of threads to use
# @param outputDir Working directory where it must save the output
#

consensusBAM <- commandArgs(trailingOnly = T)[1]
nThreads <- commandArgs(trailingOnly = T)[2]
outputDir <- commandArgs(trailingOnly = T)[3]
samtoolsDir <- commandArgs(trailingOnly = T)[4]
bedtoolsDir <- commandArgs(trailingOnly = T)[5]

# Modify the software path
if(bedtoolsDir!="none"){
    bedtoolsDir <- paste0(bedtoolsDir, "/")
}else{
    bedtoolsDir <- ""
}
if(samtoolsDir!="none"){
    samtoolsDir <- paste0(samtoolsDir, "/")
}else{
    samtoolsDir <- ""
}

# Load packages
library(here)

# Load custom functions
source(paste0(here(), "/Fastq_Functions.R"))

main <- function(){
    #Set the working directory
    setwd(outputDir)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = paste0("\nRunning the Fastq module:\n1. Filter the consensus BAM file"), 
               outDir = outputDir)
    
    # Filter the consensus BAM file
    filterBAM(bam = consensusBAM, nCores = nThreads, 
              filteredBamName = "sorted_filtered.bam",
              samtoolsDir = samtoolsDir)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "2. Create the FASTQ files", 
               outDir = outputDir)
    
    # Create the FASTQ files
    fastqCreator(bam = "sorted_filtered.bam", bedtoolsDir = bedtoolsDir)
    
}

main()