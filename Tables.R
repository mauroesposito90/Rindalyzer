## The function of the module table is to create the tables for the filtering 
## from the BAM
# @param bam Link to BAM file
# @param outputDir DIrectory for the output of the table
# @param nThreads Number of threads to use
# @param samName Name of the resulting SAM
# @param typeSam Is it a SAM derived from the mapping on consensus or genome
#

bam <- commandArgs(trailingOnly = T)[1]
outputDir <- commandArgs(trailingOnly = T)[2]
nThreads <- as.numeric(commandArgs(trailingOnly = T)[3])
samName <- commandArgs(trailingOnly = T)[4]
typeSam <- commandArgs(trailingOnly = T)[5]
samtoolsDir <- commandArgs(trailingOnly = T)[6]

# Load packages
library(here)

# Load custom functions
source(paste0(here(), "/Tables_Functions.R"))

# Modify the software path
if(samtoolsDir!="none"){
    samtoolsDir <- paste0(samtoolsDir, "/")
}else{
    samtoolsDir <- ""
}

# Constants
CHUNK <- 10000000

# Define the number of threads to use
setDTthreads(nThreads)

main <- function(){
    #Set the working directory
    setwd(outputDir)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = paste0("\nRunning the Tables module on ", typeSam, 
                                " BAM:\n1. Filtering the BAM and transforming to SAM"), 
               outDir = outputDir)
    
    # Transform BAM to SAM and filter it
    system(paste0(samtoolsDir, "samtools view -F 2308 -@ ", nThreads, " -o ", 
                  samName, " ",bam))
    
    # Compute the number of lines of the SAM file
    nLines <- nLinesSAM(sam = paste0(samName))
    
    # Check the SAM file
    if(nLines==0){
        # Update the log file
        logUpdater(modality = "Update", 
                   message = "ERROR. No reads after the filtering.", 
                   outDir = outputDir)
        
        quit(save = "no")
    }
    
    # Calculate where the SAM file must be split
    startCut <- seq(1, nLines, CHUNK)
    
    # Define the counter for the log
    counter <- 2:(length(startCut)*10+1)
    
    # Split the SAM file and work on it
    for(eachStart in 1:length(startCut)){
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], ". Open on the SAM block ", eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Open the splitted SAM
        splittedSam <- openSplitSAM(rowCut = startCut, index = eachStart, 
                                    sam = samName, SAMtype = typeSam, 
                                    nCores = nThreads)
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Compute the number of mismatches on the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Compute the number of mismatches for each read
        splittedSam[,NM:=mismCounter(NMstring = splittedSam[,NM])]
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Compute the number of H/S on the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Compute the number of H/S for each read
        splittedSam[,totHS:=HScounter(CIGAR), by = seq_len(nrow(splittedSam))]
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Compute the number of mismatches+H/S on the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Compute the total number of mismatches/H/S for each read
        splittedSam[,':='(totMism=NM+totHS,
                          NM=NULL,
                          totHS=NULL)]
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Identify the real start of reads on the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Identify the real start of each read considering H/S
        splittedSam[,realStart:=realStart(cigar = CIGAR, starting = start), 
                    by = seq_len(nrow(splittedSam))]
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Identify the real end of reads on the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Identify the real end of each read considering H/S
        splittedSam[,realEnd:=realEnd(cigar = CIGAR, starting = start), 
                    by = seq_len(nrow(splittedSam))]
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Identify the putative end of reads on the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Identify the putative end of each read without considering H/S
        splittedSam[,end:=putativEnd(cigar = CIGAR, starting = start), 
                    by = seq_len(nrow(splittedSam))]
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Identify the strand of reads on the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Identify the strand
        splittedSam[,strand:=strandIfier(flag = FLAG), 
                    by = seq_len(nrow(splittedSam))]
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Reshape the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        # Reshape the table
        splittedSam <- reShaper(tabToReshape = splittedSam, type = typeSam)
        
        # Update the log file and counter
        logUpdater(modality = "Update", 
                   message = paste0(counter[1], 
                                    ". Save the SAM block ", 
                                    eachStart), 
                   outDir = outputDir)
        counter <- counter[-1]
        
        ### Save the table ###
        writerTab(type = typeSam, tabToSave = splittedSam, nCore = nThreads)
        
    }
}

main()
