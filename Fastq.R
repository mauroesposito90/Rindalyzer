# Copyrigth (C) 2023 Mauro Esposito

# This file is part of Rindalyzer.

# Rindalyzer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Rindalyzer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Rindalyzer.  If not, see <http://www.gnu.org/licenses/>.

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