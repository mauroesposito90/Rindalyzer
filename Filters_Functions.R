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

## The set of custom functions created for the filters module

# Load packages
library(data.table)
library(seqinr)
library(stringr)

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

## The openTable function open the table consensus/genome
# @param fileTable Link to the consensus table
# @param tabType Table consensus or genome?
# @param nCores Number of threads to use
# @return tab Consensus table
#
openTable <- function(fileTable, nCores, tabType){
    tab <- fread(input = fileTable, sep = "\t", header = T, 
                 stringsAsFactors = F, nThread = nCores, 
                 colClasses = c(readName="character",realStart="numeric",
                                realEnd="numeric", totMism="numeric",
                                readSequence="character"), 
                 drop = c(2,3,4,8), key = "readName")
    return(tab)
}

## The consLen function identifies the length of the consensus sequence
# @param consRef Link to consensus reference
# @return length Length of the consensus sequence
#
consLen <- function(consRef){
    fasta <- read.fasta(file = consRef)
    length <- getLength(object = fasta)
    return(length)
}

## The fragmentIdentifierConsensus function identifies the Inside/Outside fragments
# @param tabConsensus Consensus table
# @param consensusLength Length of the consensus sequence
# @param typeFragment Type of fragment to identify (Inside or Outside?)
# @return fragments Table with identified fragments
#
fragmentIdentifierConsensus <- function(tabConsensus, consensusLength, typeFragment){
    if(typeFragment=="Inside"){
        if(nrow(tabConsensus[,if(.N==2){sum(.N)},by=readName])!=0){
            fragments <- data.table("readName"=tabConsensus[,if(.N==2){.SD}, by=readName]
                                    [realStart>0 & realEnd<=consensusLength, 
                                        if(.N==2){readName}, 
                                        by=readName]
                                    [,readName]
                                    , stringsAsFactors = F)
        }else{
            fragments <- as.data.table(matrix(nrow = 0, ncol = 1, 
                                              dimnames = list(NULL, "readName")))
        }
    }else if(typeFragment=="Outside"){
        if(nrow(tabConsensus[,if(.N==1){.SD}, by=readName])!=0){
            fragments <- data.table("readName"=tabConsensus[,if(.N==1){.SD}, by=readName]
                                    [realStart>0 & realEnd<=consensusLength, readName],
                                    stringsAsFactors = F)
        }else{
            fragments <- as.data.table(matrix(nrow = 0, ncol = 1, 
                                              dimnames = list(NULL, "readName")))
        }
    }
    return(fragments)
}

## The createRes function creates the empty table with the results of the filtering
# @return newResTable Empty table with the results to create
#
createRes <- function(){
    newResTable <- data.frame("Read"=c("Inside", "Inside"), "Mate"=c("Inside", "Outside"),
                              "ConsensusPreFilter"=c(0,0), "ConsensusPostFilter"=c(0,0),
                              "GenomePreFilter"=c(0,0), "GenomePostFilter"=c(0,0))
    return(newResTable)
}

## The updateRes function update the results of the filtering
# @param insFragments Table with the inside fragments
# @param outFragments Table with the outside fragments
# @param filteringStep Step of the filter (consPre, consPost, genPre, genPost)
# @param resTable Table with the results to update
# @return newResTable Table with the results to create or to update
#
updateRes <- function(insFragments, outFragments, filteringStep, resTable){
    if(filteringStep=="consPre"){colToUpdate <- 3
    }else if(filteringStep=="consPost"){colToUpdate <- 4
    }else if(filteringStep=="genPre"){colToUpdate <- 5
    }else if(filteringStep=="genPost"){colToUpdate <- 6}
    newResTable <- resTable
    newResTable[, colToUpdate] <- c(nrow(insFragments), nrow(outFragments))
    return(newResTable)
}

## The updateReads function update the reads information after the filtering
# @param insFragments Table with the inside fragments
# @param outFragments Table with the outside fragments
# @param filteringStep Step of the filter (ConsensusPreFilter, ConsensusPostFilter, 
#        GenomePreFilter, GenomePostFilter)
# @param readsTable Table with the reads to update
# @return newReadsTable Table with the reads to update
#
updateReads <- function(filteringStep, readsTable, insFragments, outFragments){
    readsTable[, filteringStep] <- "Absent"
    readsTable[readsTable$readName %in% insFragments$readName, filteringStep] <- "insIns"
    readsTable[readsTable$readName %in% outFragments$readName, filteringStep] <- "insOut"
    newReadsTable <- readsTable
    return(newReadsTable)
}

## The filterTable function apply the filter on the consensus/genome table
# @param fragments Table with the fragments mapped on the consensus/genome
# @param nMaxMism Maximal percentage of mismatches allowed 
# @return newTable Table with the reads filtered
#
filterTable <- function(fragments, nMaxMism){
    newTable <- fragments
    newTable[,readLength:=str_count(readSequence)]
    newTable[,maxMism:=floor(readLength*nMaxMism/100)]
    newTable <- newTable[!(newTable$readName %in% newTable[totMism>maxMism, readName]),]
    newTable[,":="(readLength=NULL, maxMism=NULL)]
    return(newTable)
}

## The fragmentIdentifierGenome function identifies the Inside/Outside fragments
# @param tabGenome Genome table
# @param typeFragment Type of fragment to identify (insIns or insOut?)
# @param readsInfo Information about each reads
# @return fragments Table with identified fragments
#
fragmentIdentifierGenome <- function(tabGenome, typeFragment, readsInfo){
    if(typeFragment=="insIns"){
        insideConsensus <- unique(readsInfo[ConsensusPostFilter==typeFragment, readName])
        if(length(insideConsensus)!=0){
            fragments <- data.table("readName"=insideConsensus[insideConsensus %in% tabGenome$readName]
                                    , stringsAsFactors = F)
        }else{
            fragments <- as.data.table(matrix(nrow = 0, ncol = 1, 
                                              dimnames = list(NULL, "readName")))
        }
    }else if(typeFragment=="insOut"){
        outsideConsensus <- unique(readsInfo[ConsensusPostFilter==typeFragment, readName])
        if(length(outsideConsensus)!=0){
            fragments <- data.table("readName"=outsideConsensus[outsideConsensus %in% tabGenome$readName]
                                    , stringsAsFactors = F)
        }else{
            fragments <- as.data.table(matrix(nrow = 0, ncol = 1, 
                                              dimnames = list(NULL, "readName")))
        }
    }
    return(fragments)
}

## The bedCreator function creates the BED file of the genomic reads to intersect
# @param tabGenome Genome table 
# @param readsInfo Information about each reads
# @param nCores Number of threads to use
#
bedCreator <- function(tabGenome, readsInfo, nCores){
    readsInfo <- readsInfo[GenomePreFilter!="Absent", c("readName", "GenomePreFilter")]
    readsInfo <- readsInfo[!duplicated(readsInfo),]
    tabGenome <- merge.data.table(x = tabGenome, y = readsInfo, by = "readName")
    tabGenome <- tabGenome[,c("chromosome", "realStart", "realEnd", "readName", "GenomePreFilter")]
    fwrite(x = tabGenome, file = "fragments.bed", sep = "\t", row.names = F, col.names = F, 
           nThread = nCores)
}

## The nLines function calculate the number of lines of a file
# @param file Link to the file
# @return lines Number of lines
#
nLines <- function(file){
    wcLines <- system(paste0("wc -l ", file), intern = T)
    lines <- as.numeric(strsplit(wcLines, " ")[[c(1,1)]])
    return(lines)
}

## The openIntersections function open the file with the reads-BED intersections
# @param file Link to the file of intersections
# @param nCores Number of threads to use
# @return table Table with the intersections
#
openIntersections <- function(file, nCores){
    table <- fread(file = file, sep = "\t", stringsAsFactors = F, header = F, 
                   nThread = nCores, col.names = c("chromosome", "start", "end", 
                                                   "readName", "type", "chromosome2", 
                                                   "start2", "end2", "nOverlap"))
}

## The filterIntersections function filter out the intersections based on the 
## maximum number of mismatches allowed
# @param intersections Table of intersections
# @param percentageMismatches Percentage of mismatches allowed
# @return newIntersections Table with the filtered intersections
#
filterIntersections <- function(intersections, percentageMismatches){
    newIntersections <- intersections
    newIntersections[,length:=end-start]
    newIntersections[,minOverlap:=length-floor(length/100*percentageMismatches)]
    newIntersections <- newIntersections[nOverlap>=minOverlap,]
    return(newIntersections)
}
## The fragmentIdentifierBed function identifies the Inside/Outside fragments 
## after the BED filter
# @param intersections Table of intersections
# @param typeFragment Type of fragment to identify (insIns or insOut?)
# @return fragments Table with identified fragments
#
fragmentIdentifierBed <- function(intersections, typeFragment){
    fragments <- unique(intersections[type==typeFragment, readName])
    if(length(fragments)!=0){
        fragments <- data.table("readName"=fragments, stringsAsFactors = F)
    }else{
        fragments <- as.data.table(matrix(nrow = 0, ncol = 1, 
                                          dimnames = list(NULL, "readName")))
    }
    return(fragments)
}
