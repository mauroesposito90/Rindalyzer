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
library(argparser)

## The testNoInput function stop the code if there are no input
#
testNoInput <- function(){
    commandLine <- commandArgs(trailingOnly = T)
    nArguments <- length(grepl(pattern = "[:graph:]", x = commandLine))
    if(nArguments==0){
        print("No arguments provided. Please check ./Rindalyzer --help for further information")
        quit()
    }
}

## The parserCreator function create parser object
# @return parser Parser object
#
parserCreator <- function(){
    parser <- arg_parser(name = "Rindalyzer", 
                         description = "Rindalyzer tool allows the to highlight active transcription of TEs", 
                         hide.opts = T)
    parser <- add_argument(parser = parser, arg = "--FASTQ_1", 
                           help = "Path to FASTQ r1 file", 
                           type = "character", short = "-1", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--FASTQ_2",
                           help = "Path to FASTQ r2 file",
                           type = "character", short = "-2", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--CONSENSUS",
                           help = "Path to FASTA file of consensus sequence",
                           type = "character", short = "-c", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--MISMATCHES",
                           help = "Maximum percentage of allowed mismatches",
                           type = "numeric", short = "-m", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--GENOME",
                           help = "Path to FASTA file of the reference genome",
                           type = "character", short = "-g", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--OUTPUT",
                           help = "Path to directory for the output",
                           type = "character", short = "-o", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--NAME",
                           help = "Name of the analysis",
                           type = "character", short = "-n", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--THREADS",
                           help = "Number of threads",
                           type = "numeric", short = "-t", default = 1, nargs = 1)
    parser <- add_argument(parser = parser, arg = "--BED",
                           help = "Path to BED file of genomic regions",
                           type = "character", short = "-b", default = "none", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--KEEP_TEMPORARY",
                           help = "TRUE for maintaining temporary files",
                           type = "logic", short = "-k", default = FALSE, nargs = 1)
    parser <- add_argument(parser = parser, arg = "--SAMTOOLS",
                           help = "Path to directory containing the Samtools software",
                           type = "character", short = "-s", default = "none", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--BEDTOOLS",
                           help = "Path to directory containing the Bedtools software",
                           type = "character", short = "-B", default = "none", nargs = 1)
    parser <- add_argument(parser = parser, arg = "--BWA",
                           help = "Path to directory containing the BWA software",
                           type = "character", short = "-w", default = "none", nargs = 1)
    parser <- parse_args(parser = parser)
    return(parser)
}

## The normalizator function create normalized paths in the parser
# @param parserTemp A parser
# @return parserTemp A parser with the normalized paths
#
normalizator <- function(parserTemp){
    for(arg in 3:length(parserTemp)){
        if(grepl(pattern = "[/]", x = unlist(parserTemp[arg]))){
            parserTemp[arg] <- normalizePath(unlist(parserTemp[arg]))
        }
    }
    return(parserTemp)
}

## The logUpdater function create/update the log file
# @param modality Create or Update the log file
# @param message Message to write in the log file
# @param outDir Directory for the output
#
logUpdater <- function(modality, message, outDir){
    if(modality=="Create"){
        file.create("Rindalyzer.log")
        appendOption <- FALSE
    }else if(modality=="Update"){
        appendOption <- TRUE
    }
    write(x = message, file = paste0(outDir, "/Rindalyzer.log"), 
          append = appendOption)
}

## The tempDelete function delete the temporary files
# @param parser Parser of the code
#
tempDelete <- function(parser){
    if(parser$KEEP_TEMPORARY==FALSE){
        unlink(c(paste0(parser$OUTPUT, "/tempSam"),
                 paste0(parser$OUTPUT, "/rindGenome*"),
                 paste0(parser$OUTPUT, "/rindConsensus*"),
                 paste0(parser$OUTPUT, "/*.bam"),
                 paste0(parser$OUTPUT, "/*.fastq")))
        if(parser$BED!="none"){
            unlink(c(paste0(parser$OUTPUT, "/intersected.txt"),
                     paste0(parser$OUTPUT, "/fragments.bed")))
        }
    }
}
