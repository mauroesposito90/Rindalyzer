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

## The function of the module mapping is to map FASTQ files on the consensus 
## and on the genome
# @param fastq_r1 Link to FASTQ files
# @param fastq_r2 Link to FASTQ files
# @param reference Link to Reference file
# @param indexName Name of the BWA index
# @param bamName Name of the resulting BAM
# @param nThreads Number of threads to use
# @param outputDir Directory for the index and BAM file
# @param type Consensus or genome?
#

fastq_r1 <- commandArgs(trailingOnly = T)[1]
fastq_r2 <- commandArgs(trailingOnly = T)[2]
reference <- commandArgs(trailingOnly = T)[3]
indexName <- commandArgs(trailingOnly = T)[4]
bamName <- commandArgs(trailingOnly = T)[5]
nThreads <- commandArgs(trailingOnly = T)[6]
outputDir <- commandArgs(trailingOnly = T)[7]
type <- commandArgs(trailingOnly = T)[8]
bwaDir <- commandArgs(trailingOnly = T)[9]

# Load packages
library(here)

# Load custom functions
source(paste0(here(), "/Mapping_Functions.R"))

# Modify the software path
if(bwaDir!="none"){
    bwaDir <- paste0(bwaDir, "/")
}else{
    bwaDir <- ""
}

main <- function(){
    # Set the working directory
    setwd(outputDir)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = paste0("\nRunning the Mapping module (", type, 
                                "):\n1. Check the FASTQ format"), 
               outDir = outputDir)
    
    # Check the FASTQ format of the files
    validityFastq_r1 <- fastqValidator(fastq = fastq_r1)
    validityFastq_r2 <- fastqValidator(fastq = fastq_r2)
    if(validityFastq_r1==F | validityFastq_r2==F){
        # Update the log file
        logUpdater(modality = "Update", 
                   message = "ERROR. Invalid FASTQ format.", 
                   outDir = outputDir)
        
        quit(save = "no")
    }
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = paste0("2. Check the FASTA format of ", type, " file"), 
               outDir = outputDir)
    
    # Check the FASTA format of the Reference
    validityReference <- fastaValidator(fasta = reference)
    if(validityReference==F){
        # Update the log file
        logUpdater(modality = "Update", 
                   message = paste0("ERROR. Invalid ", type, " FASTA format."), 
                   outDir = outputDir)
        
        quit(save = "no")
    }
    
    # Define the BWA algorithm for the indexing step
    BWAalgorithm <- indexAlgorithm(reference = reference)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = paste0("3. Indexing the ", type, " file"), 
               outDir = outputDir)
    
    # Run BWA index
    system(paste0(bwaDir, "bwa index -p ", indexName, " -a ", BWAalgorithm, " ",reference))
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = paste0("4. Map the FASTQ on the ", type, " FASTA"), 
               outDir = outputDir)
    
    # Run mapping step
    system(paste0(bwaDir, "bwa mem -t ", nThreads, " ", indexName, " ", fastq_r1, " ",
                 fastq_r2, " > ", bamName, ".bam"))
}

main()