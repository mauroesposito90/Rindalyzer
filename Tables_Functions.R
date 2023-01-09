## The set of custom functions created for the tables module

# Load packages
library(data.table)
library(Rsamtools)
library(stringr)
library(GenomicAlignments)

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

## The openSplitSAM function opens some lines of a file starting from the beginning
# @param rowCut Number of the lines for starting the cutting
# @param index Index for choosing the rowCut
# @param sam Link to the SAM file to open
# @param nCores Number of threads
# @param SAMtype Type of SAM (consensus or genomic)
# @return splittedSam Table with the splitted SAM
#
openSplitSAM=function(rowCut, index, sam, SAMtype, nCores){
    rowsToSkip <- ifelse(test = rowCut[index]==1, yes = 0, no = rowCut[index]-1)
    rowsToOpen <- ifelse(test = index!=length(rowCut), no = Inf,
                        yes = rowCut[index+1]-rowCut[index])
    if(SAMtype=="consensus"){
        colsToSelect <- c("1,4,6,2,10")
        colsClasses <- c("character", rep("numeric",2), rep("character",2))
        colsNames <- c("readName", "FLAG","start", "CIGAR", "readSequence")
    }else if(SAMtype=="genome"){
        colsToSelect <- c("1,4,6,2,10,3")
        colsClasses <- c(rep(c("character","numeric"),2), rep("character",2))
        colsNames <- c("readName", "FLAG", "chromosome","start", "CIGAR", "readSequence")
    }
    splittedSam <- fread(cmd = paste0("cut -f ", colsToSelect, " ", sam), sep = "\t", 
                      nThread = nCores, nrows = rowsToOpen, skip = rowsToSkip, 
                      stringsAsFactors = F, colClasses = colsClasses, header = F, 
                      col.names = colsNames)
    NMinfo <- fread(cmd = paste0("grep -o NM:i:[0-9]* ", sam), 
                       nThread = nCores, header = F, nrows = rowsToOpen,
                       skip = rowsToSkip, stringsAsFactors = F, 
                       colClasses = c("character"), col.names = "NM")
    splittedSam <- cbind(splittedSam, NMinfo)
    return(splittedSam)
}

## The strandInfier function identify the strand of the mapped reads
# @param flag FLAG of the read
# @return strand Strand of the mapped read
#
strandIfier <- function(flag){
    strand <- as.numeric(as.vector(bamFlagAsBitMatrix(as.integer(flag))[1,5]))
    return(strand)
}

## The mismCounter function sum the number of mismatches based on the NM info
# @param NMstring NM information of the mapped read
# @return nMism Number of mismatches
#
mismCounter <- function(NMstring){
    nMism <- as.numeric(gsub(x = NMstring, pattern = "NM:i:", replacement = ""))
    return(nMism)
    }

## The HScounter function sum the number of hard/soft clipped nucleotides
# @param cigar CIGAR information of the mapped read
# @return totalHS Number of H+S
#
HScounter <- function(cigar){
    letterS <- unlist(str_extract_all(string = cigar, 
                                      pattern = "[0-9]+S"))
    totalS <- sum(as.numeric(gsub(x = letterS, 
                                  pattern = "[a-zA-Z]", 
                                  replacement = "")))
    letterH <- unlist(str_extract_all(string = cigar, pattern = "[0-9]+H"))
    totalH <- sum(as.numeric(gsub(x = letterH, 
                                  pattern = "[a-zA-Z]", 
                                  replacement = "")))
    totalHS <- totalS + totalH
    return(totalHS)
}

## The realStart function identifies the coordinate of the real start considering 
## the initial H/S
# @param cigar CIGAR information of the mapped read
# @param starting Start of the mapped read
# @return trueStart Real start of the read
#
realStart <- function(cigar, starting){
    softNT <- grep(x = cigar, pattern = "^[0-9]+S")
    if(length(softNT)!=0){
        nSoft <- as.numeric(gsub(x = str_extract(string = cigar, pattern = "^[0-9]+S"), 
                                 pattern = "S", 
                                 replacement = ""))
    }else{
        nSoft <- 0
    }
    hardNT <- grep(x = cigar, pattern = "^[0-9]+H")
    if(length(hardNT)!=0){
        nHard <- as.numeric(gsub(x = str_extract(string = cigar, pattern = "^[0-9]+H"), 
                                 pattern = "H", 
                                 replacement = ""))
    }else{
        nHard <- 0
    }
    trueStart <- starting-nSoft-nHard
    return(trueStart)
}

## The realEnd function identifies the coordinate of the real end considering 
## the final H/S
# @param cigar CIGAR information of the mapped read
# @param starting Start of the mapped read
# @return trueEnd Real end of the read
#
realEnd <- function(cigar, starting){
    softNT <- grep(x = cigar, pattern = "[0-9]+S$")
    if(length(softNT)!=0){
        nSoft <- as.numeric(gsub(x = str_extract(string = cigar, pattern = "[0-9]+S$"), 
                                 pattern = "S", 
                                 replacement = ""))
    }else{
        nSoft <- 0
    }
    hardNT <- grep(x = cigar, pattern = "[0-9]+H$")
    if(length(hardNT)!=0){
        nHard <- as.numeric(gsub(x = str_extract(string = cigar, pattern = "[0-9]+H$"), 
                                 pattern = "H", 
                                 replacement = ""))
    }else{
        nHard <- 0
    }
    inner <- cigarWidthAlongReferenceSpace(cigar)
    trueEnd <- inner+nSoft+nHard+starting-1
    return(trueEnd)
}

## The putativEnd function identifies the coordinate of the end without considering 
## the final H/S
# @param cigar CIGAR information of the mapped read
# @param starting Start of the mapped read
# @return falsEnd Putative end of the read
#
putativEnd <- function(cigar, starting){
    inner <- cigarWidthAlongReferenceSpace(cigar)
    falsEnd <- inner+starting-1
    return(falsEnd)
}

## The nLinesSAM function calculate the number of lines of a SAM file
# @param sam Link to the SAM file
# @return lines Number of lines
#
nLinesSAM <- function(sam){
    wcLines <- system(paste0("wc -l ", sam), intern = T)
    lines <- as.numeric(strsplit(wcLines, " ")[[c(1,1)]])
    return(lines)
}

## The writerTab function save the final table
# @param type The tables mapping is running on consensus or genomic?
# @param tabToSave Table to save
# @param nCore Number of threads
#
writerTab <- function(type, tabToSave, nCore){
    fileToSave <- paste0("./table_", type, ".txt")
    appendOption <- ifelse(test = file.exists(fileToSave), yes = T, no = F)
    colnamesOption <- ifelse(test = file.exists(fileToSave), yes = F, no = T)
    fwrite(x = tabToSave, file = fileToSave, append = appendOption, quote = F, 
           sep = "\t", nThread = nCore, col.names = colnamesOption)
}

## The reShaper function reshape the final table
# @param tabToReshape Table to reshape
# @param type The tables mapping is running on consensus or genomic?
# @return shapedTab Table ordered
#
reShaper <- function(tabToReshape, type){
    tabToReshape[,FLAG:=NULL]
    if(type=="consensus"){
        cols <- c("readName", "start", "end", "CIGAR", "realStart", 
                  "realEnd", "totMism", "strand", "readSequence")
    }else if(type=="genome"){
        cols <- c("readName", "start", "end", "CIGAR", "realStart", 
                  "realEnd", "totMism", "strand", "readSequence", "chromosome")
    }
    setcolorder(x = tabToReshape, neworder = cols)
    return(tabToReshape)
}
