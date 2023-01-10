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