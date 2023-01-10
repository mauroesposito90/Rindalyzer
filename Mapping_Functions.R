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

## The set of custom functions created for the mapping module

# Load packages
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

## The openFile function open some lines of a file starting from the beginning
# @param file Link to the file to open
# @param nLines Number of lines to open
# @return lines Lines read
#
openFile <- function(file, nLines){
    formatNotCompressed <- c(".fa", "fasta", "fastq", ".fq")
    if(str_ends(string = file, pattern = ".gz")){
        lines <- system(command = paste0("zcat ", file, " 2>/dev/null | head -n ", 
                                         nLines), 
                        intern = T)
    }else if (sum(str_ends(string = file, pattern = formatNotCompressed))!=0){
        lines <- system(command = paste0("head -n ", nLines, " ",file), 
                        intern = T)
    }else{
        print("File not supported")
        quit(save = "no")
    }
    return(lines)
}

## The fastqValidator function validates the format of the FASTQ file to map
#  by checking the firsts 8 lines 
# @param fastq Link to the FASTQ file to validate
# @return valid True or False based on the controls of the FASTQ file
#
fastqValidator <- function(fastq){
    firstLines <- openFile(fastq, 8)
    
    firstLineCheck <- str_starts(string = firstLines[c(1, 5)], pattern = "@")
    thirdLineCheck <- str_starts(string = firstLines[c(3, 7)], pattern = "+")
    secondFourthLinesCheck <- c(length(firstLines[2]) == length(firstLines[4]),
                                length(firstLines[6]) == length(firstLines[8]))
    valid <- ifelse(test = all(firstLineCheck, thirdLineCheck, secondFourthLinesCheck),
                    yes = valid <- T, 
                    no = valid <- F)
    return(valid)
}

## The referenceValidator function validates the FASTA format of the reference file
#  by checking the firsts 2 lines 
# @param fasta Link to the FASTA file to validate
# @return valid True or False based on the controls of the FASTA file
#
fastaValidator <- function(fasta){
    firstLines <- openFile(file = fasta, nLines = 2)
    firstLineCheck <- str_starts(string = firstLines[1], pattern = ">")
    secondLineCheck <- str_starts(string = firstLines[2], pattern = c("[A-Z]", 
                                                                      "[a-z]"))
    valid <- ifelse(test = sum(firstLineCheck, secondLineCheck)==2,
                    yes = T, 
                    no = F)
    return(valid)
}

## The indexAlgorithm function decide the type of indexing for BWA based on 
## the size of the reference
# @param reference Link to the reference file
# @return algorithm Name of the algorithm to use
#
indexAlgorithm <- function(reference){
    fileSize <- file.info(reference)$size
    algorithm <- ifelse(test = fileSize<2e+09, yes = "is", no = "bwtsw")
    return(algorithm)
}
