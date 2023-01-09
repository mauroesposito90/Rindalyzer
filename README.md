# Rindalyzer
Tool for the autonomous transcription of transposable elements

# Overview
Rindalyzer is a tool under development that can highlight the autonomous transcription of transposable elements when these result overexpressed. Check the Esposito et al., Biomedicines, 2022 (https://doi.org/10.3390/biomedicines10123279) pubblication for more details.

# Prerequisites
To run Rindalyzer you need of the following prerequisites:
* R packages: data.table, seqinr, here, stringr, ggplot2, tidyr, dplyr, magrittr, argparser, Rsamtools and GenomicAlignments
* Software: Samtools, Bedtools, BWA

# How run Rindalyzer
The tool can be used from the command-line by running: “./Rindalyzer.R [options]”. 

The options are:
* -1 / --FASTQ_1, Path to FASTQ r1 file
* -2 / --FASTQ_2, Path to FASTQ r2 file
* -c / --CONSENSUS, Path to FASTA file of consensus sequence
* -m / --MISMATCHES, Maximum percentage of allowed mismatches
* -g / --GENOME, Path to FASTA file of the reference genome
* -o / --OUTPUT, Path to directory for the output
* -n / --NAME, Name of the analysis
* -t / --THREADS, Number of threads [default: 1]
* -b / --BED, Path to BED file of genomic regions [default: none]
* -k / --KEEP_TEMPORARY, TRUE for maintaining temporary files [default: FALSE]
* -s / --SAMTOOLS, Path to directory containing the Samtools software [default: none]
* -B / --BEDTOOLS, Path to directory containing the Bedtools software [default: none]
* -w / --BWA, Path to directory containing the BWA software [default: none]
 
Example usage: 
./Rindalyzer.R --FASTQ_1 /path/to/fastq1.fastq --FASTQ_2 /path/to/fastq2.fastq --CONSENSUS /path/to/consensus.fa --BED /path/to/file.bed --MISMATCHES 10 --GENOME /path/to/genome.fasta --OUTPUT /path/to/output/ --NAME test [ … ]

The output is reported in the file “rindReport.html”.
