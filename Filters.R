## The function of the module filters is to count the reads based on the filters
# @param consensusTable Link to consensus table
# @param genomeTable Link to genomic table
# @param consensusReference Link to reference consensus sequence
# @param mismatchPercentage Maximal percentage of mismatches allowed
# @param nThreads Number of threads to use
# @param outputDir Directory for the index and BAM file
# @param bed Link to the BED file for the filter
#

consensusTable <- commandArgs(trailingOnly = T)[1]
genomeTable <- commandArgs(trailingOnly = T)[2]
consensusReference <- commandArgs(trailingOnly = T)[3]
mismatchPercentage <- as.numeric(commandArgs(trailingOnly = T)[4])
nThreads <- as.numeric(commandArgs(trailingOnly = T)[5])
outputDir <- commandArgs(trailingOnly = T)[6]
bed <- commandArgs(trailingOnly = T)[7]
bedtoolsDir <- commandArgs(trailingOnly = T)[8]

# Load packages
library(here)

# Load custom functions
source(paste0(here(), "/Filters_Functions.R"))

# Modify the software path
if(bedtoolsDir!="none"){
    bedtoolsDir <- paste0(bedtoolsDir, "/")
}else{
    bedtoolsDir <- ""
}

main <- function(){
    # Set the working directory
    setwd(outputDir)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "\nRunning the Filters module:\n1. Open the consensus table", 
               outDir = outputDir)
    
    # Open the consensus table
    tableConsensus <- openTable(fileTable = consensusTable, nCores = nThreads)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "2. Identify the Inside/Outside fragments", 
               outDir = outputDir)
    
    # Identify the length of the consensus sequence
    lengthConsensus <- consLen(consRef = consensusReference)
    
    # Identify the Inside/Outside fragments
    insideFragments <- fragmentIdentifierConsensus(tabConsensus = tableConsensus, 
                                                   consensusLength = lengthConsensus, 
                                                   typeFragment = "Inside")
    outsideFragments <- fragmentIdentifierConsensus(tabConsensus = tableConsensus, 
                                                    consensusLength = lengthConsensus, 
                                                    typeFragment = "Outside")
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "3. Create the consensusPreFilter results", 
               outDir = outputDir)
    
    # Create the table with the consensus pre-filter results
    results <- createRes()
    results <- updateRes(insFragments = insideFragments, outFragments = outsideFragments, 
                         filteringStep = "consPre", resTable = results)
    
    # Update the table with the information about each single read
    readsTable <- updateReads(filteringStep = "ConsensusPreFilter", readsTable = tableConsensus, 
                              insFragments = insideFragments, outFragments = outsideFragments)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "4. Apply the mismatch filter", 
               outDir = outputDir)
    
    # Apply the mismatch filter on the consensus mapping
    tableConsensus <- filterTable(fragments = tableConsensus, 
                                  nMaxMism = mismatchPercentage)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "5. Identify the Inside/Outside fragments", 
               outDir = outputDir)
    
    # Identify the Inside/Outside fragments
    insideFragments <- fragmentIdentifierConsensus(tabConsensus = tableConsensus, 
                                                   consensusLength = lengthConsensus, 
                                                   typeFragment = "Inside")
    outsideFragments <- fragmentIdentifierConsensus(tabConsensus = tableConsensus, 
                                                    consensusLength = lengthConsensus, 
                                                    typeFragment = "Outside")
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "6. Create the consensusPostFilter results", 
               outDir = outputDir)
    
    # Update the table with the consensus post-filter results
    results <- updateRes(insFragments = insideFragments, outFragments = outsideFragments, 
                         filteringStep = "consPost", resTable = results)
    
    # Update the table with the information about each single read
    readsTable <- updateReads(filteringStep = "ConsensusPostFilter", readsTable = readsTable,
                              insFragments = insideFragments, outFragments = outsideFragments)
    rm(tableConsensus)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "7. Open the genome table", 
               outDir = outputDir)
    
    # Open the genome table
    tableGenome <- openTable(fileTable = genomeTable, nCores = nThreads)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "8. Filter the table based on the consensus mapping", 
               outDir = outputDir)
    
    # Discard the reads that have not passed the consensus filters
    tableGenome <- tableGenome[readName %in% readsTable[ConsensusPostFilter!="Absent", readName]]
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "9. Identify the Inside/Outside fragments", 
               outDir = outputDir)
    
    # Identify the Inside/Outside fragments
    insideFragments <- fragmentIdentifierGenome(tabGenome = tableGenome, 
                                                typeFragment = "insIns", 
                                                readsInfo = readsTable)
    outsideFragments <- fragmentIdentifierGenome(tabGenome = tableGenome, 
                                                 typeFragment = "insOut", 
                                                 readsInfo = readsTable)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "10. Create the genomePreFilter results", 
               outDir = outputDir)
    
    # Update the table with the genome pre-filter results
    results <- updateRes(insFragments = insideFragments, outFragments = outsideFragments, 
                         filteringStep = "genPre", resTable = results)
    
    # Update the table with the information about each single read
    readsTable <- updateReads(filteringStep = "GenomePreFilter", readsTable = readsTable,
                              insFragments = insideFragments, outFragments = outsideFragments)
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = "11. Apply the mismatch filter", 
               outDir = outputDir)
    
    # Apply the mismatch filter on the genome mapping
    tableGenome <- filterTable(fragments = tableGenome, nMaxMism = mismatchPercentage)
    
    # Apply the BED filter on the genome mapping
    if(bed!="none"){
        # Update the log file
        logUpdater(modality = "Update", 
                   message = "12. Create the BED file of reads", 
                   outDir = outputDir)
        
        # Create the BED file of the reads
        bedCreator(tabGenome = tableGenome, readsInfo = readsTable, nCores = nThreads)
        
        # Update the log file
        logUpdater(modality = "Update", 
                   message = "13. Intersect the reads with the provided BED file", 
                   outDir = outputDir)
        
        # Intersect the reads to the provided BED file
        system(command = paste0(bedtoolsDir, "bedtools intersect -a fragments.bed -b ",
                                bed, " -wo > intersected.txt"))
        
        # Test if there are intersections
        nIntersections <- nLines(file = "intersected.txt")
        if(nIntersections!=0){
            # Update the log file
            logUpdater(modality = "Update", 
                       message = "14. Open the intersections file", 
                       outDir = outputDir)
            
            # Open the file of the intersections
            intersections <- openIntersections(file = "intersected.txt", nCores = nThreads)
            
            # Update the log file
            logUpdater(modality = "Update", 
                       message = "15. Apply the filter", 
                       outDir = outputDir)
            
            # Filter out the intersections less than maximum number of allowed mismatches
            intersections <- filterIntersections(intersections = intersections, 
                                                 percentageMismatches = mismatchPercentage)
            
            # Update the log file
            logUpdater(modality = "Update", 
                       message = "16. Identify the Inside/Outside fragments", 
                       outDir = outputDir)
            
            # Identify the Inside/Outside fragments that overcome the BED filter
            insideFragments <- fragmentIdentifierBed(intersections = intersections, 
                                                     typeFragment = "insIns")
            
            outsideFragments <- fragmentIdentifierBed(intersections = intersections, 
                                                     typeFragment = "insOut")
            
            # Create a counter
            counter <- 17
        }else{
            # Update the log file
            logUpdater(modality = "Update", 
                       message = "ERROR. No intersections with the provided BED file were found.", 
                       outDir = outputDir)
        }
    }else{
        # Update the log file
        logUpdater(modality = "Update", 
                   message = "12. Identify the Inside/Outside fragments", 
                   outDir = outputDir)
        
        # Identify the Inside/Outside fragments
        insideFragments <- fragmentIdentifierGenome(tabGenome = tableGenome, 
                                                    typeFragment = "insIns", 
                                                    readsInfo = readsTable)
        outsideFragments <- fragmentIdentifierGenome(tabGenome = tableGenome, 
                                                     typeFragment = "insOut", 
                                                     readsInfo = readsTable)
        
        # Create a counter
        counter <- 13
    }
    
    # Update the log file
    logUpdater(modality = "Update", 
               message = paste0(counter, ". Create the genomePostFilter results"), 
               outDir = outputDir)
    
    # Update the table with the genome pre-filter results
    results <- updateRes(insFragments = insideFragments, outFragments = outsideFragments, 
                         filteringStep = "genPost", resTable = results)
    
    # Update the table with the information about each single read
    readsTable <- updateReads(filteringStep = "GenomePostFilter", readsTable = readsTable,
                              insFragments = insideFragments, outFragments = outsideFragments)
 
    # Save the results
    fwrite(x = results, file = "rindResults.txt", quote = F, sep = "\t", 
           row.names = F, col.names = T, nThread = nThreads)
    fwrite(x = readsTable, file = "rindFragments.txt", quote = F, sep = "\t", 
           row.names = F, col.names = T, nThread = nThreads)
}

main()
