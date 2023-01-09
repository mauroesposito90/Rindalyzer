## The set of custom functions created for the fastq module

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

## The filterBAM function filter the consensus BAM file
# @param bam Link to the consensus BAM file
# @param nCores Number of cores to use
# @param filteredBamName Name of the resulting filtered BAM file
# @param samtoolsDir Directory of the Samtools software
#
filterBAM <- function(bam, nCores, filteredBamName, samtoolsDir){
    system(command = paste0(samtoolsDir, "samtools view -F 2304 -h -G 12 -@ ", nCores, 
                            " -o filtered.sam ", bam))
    system(command = paste0(samtoolsDir, "samtools sort -n -O BAM -@ ", nCores,
                            " -o ", filteredBamName," filtered.sam"))
    unlink("filtered.sam")
}

## The fastqCreator function creates the FASTQ files from the filtered consensus 
## BAM file
# @param bam Link to the filtered and sorted consensus BAM file
# @param bedtoolsDir Directory of the Bedtools software
#
fastqCreator <- function(bam, bedtoolsDir){
    system(command = paste0(bedtoolsDir, "bedtools bamtofastq -i ", bam, 
                            " -fq R1.fastq -fq2 R2.fastq"))
}