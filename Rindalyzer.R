#!/usr/bin/env Rscript

## The Rindalyzer tool analyzes the autonomous transcription of TEs
# @param FASTQ_1 Link to FASTQ files
# @param FASTQ_2 Link to FASTQ files
# @param CONSENSUS Link to Reference consensus file
# @param GENOME Link to Reference genome file
# @param BED Link to the BED file for the intersections
# @param MISMATCHES Maximum percentage of mismatches allowed
# @param THREADS Number of threads to use
# @param OUTPUT Directory for the index and BAM file
# @param NAME Name for the output directory
# @param KEEP_TEMPORARY TRUE for maintaining temporary files
#

# Load packages
suppressMessages(library(here))
library(argparser)

# Load custom functions
source(paste0(here(), "/Rindalyzer_Functions.R"))

main <- function(){
    # Test no input
    testNoInput()
    
    # Set the working directory
    setwd(here())
    
    # Create the parser
    parser <- parserCreator()

    # Normalize the paths in the parser
    parser <- normalizator(parser = parser)
    
    # Create the output directory
    newOutDir <- paste0(parser$OUTPUT, "/Rindalyzer_", parser$NAME)
    dir.create(newOutDir)
    parser$OUTPUT <- newOutDir

    # Create the log file
    logUpdater(modality = "Create",
               message = paste0("START Rindalyzer (", Sys.time(), ") with arguments:\n",
                                "- FASTQ_1, ", parser$FASTQ_1,
                                "\n- FASTQ_2, ", parser$FASTQ_2,
                                "\n- Consensus, ", parser$CONSENSUS,
                                "\n- Genome, ", parser$GENOME,
                                "\n- Percentage mismatches, ", parser$MISMATCHES,
                                "\n- Threads, ", parser$THREADS,
                                "\n- Analysis name, ", parser$NAME,
                                "\n- BED file, ", parser$BED,
                                "\n- Temporary files, ", parser$KEEP_TEMPORARY,
                                "\n- Output directory, ", parser$OUTPUT),
               outDir = parser$OUTPUT)

    # Map the FASTQ on the consensus
    system(command = paste("Rscript ./Mapping.R", parser$FASTQ_1, parser$FASTQ_2,
                           parser$CONSENSUS, "rindConsensusIndex rindConsensusBam",
                           parser$THREADS, parser$OUTPUT, "consensus ", parser$BWA))

    # Create the table of the consensus mapping
    system(command = paste0("Rscript ./Tables.R ", parser$OUTPUT, "/rindConsensusBam.bam ",
                            parser$OUTPUT, " ", parser$THREADS, " tempSam consensus ",
                            parser$SAMTOOLS))

    # Create the FASTQ files to map on the genome
    system(command = paste0("Rscript ./Fastq.R ", parser$OUTPUT, "/rindConsensusBam.bam ",
                            parser$THREADS, " ", parser$OUTPUT, " ", parser$SAMTOOLS,
                            " ", parser$BEDTOOLS))

    # Map the FASTQ on the genome
    system(command = paste0("Rscript ./Mapping.R ", parser$OUTPUT, "/R1.fastq ",
                            parser$OUTPUT, "/R2.fastq ", parser$GENOME,
                            " rindGenomeIndex rindGenomeBam ", parser$THREADS, " ",
                            parser$OUTPUT, " genome ", parser$BWA))

    # Create the table of the genome mapping
    system(command = paste0("Rscript ./Tables.R ", parser$OUTPUT, "/rindGenomeBam.bam ",
                            parser$OUTPUT, " ", parser$THREADS, " tempSam genome ",
                            parser$SAMTOOLS))

    # Create the table of the genome mapping
    system(command = paste0("Rscript ./Filters.R ", parser$OUTPUT, "/table_consensus.txt ",
                            parser$OUTPUT, "/table_genome.txt ", parser$CONSENSUS, " ",
                            parser$MISMATCHES, " ", parser$THREADS, " ", parser$OUTPUT,
                            " ", parser$BED, " ", parser$BEDTOOLS))

    # Create the final results
    system(command = paste0("Rscript -e \"rmarkdown::render('Report.Rmd', ",
                            "output_file='", parser$OUTPUT, "/rindReport.html',",
                            "params=list(",
                            "results='", parser$OUTPUT, "/rindResults.txt',",
                            "fragments='", parser$OUTPUT, "/rindFragments.txt',",
                            "fastq_r1='", parser$FASTQ_1, "',",
                            "fastq_r2='", parser$FASTQ_2, "',",
                            "reference='", parser$CONSENSUS, "',",
                            "bed='", parser$BED, "',",
                            "mismatchPercentage='", parser$MISMATCHES, "',",
                            "genome='", parser$GENOME, "',",
                            "reportDir='", parser$OUTPUT, "'))\""
                            )
          )

    # Discard temporary files
    tempDelete(parser = parser)
}

main()
