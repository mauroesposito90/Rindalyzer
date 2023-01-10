# Rindalyzer
This tool highlights the autonomous transcription of transposable elements (TE) when these result overexpressed. It is in development starting from a pipeline used for analyzing LINE1 elements in the human genome (Esposito et al., Biomedicines 2022, 10(12), 3279; https://doi.org/10.3390/biomedicines10123279) but the same principles can be potentially applied for the analysis of other transposable elements families and other organisms. 

# Overview
With this tool, Illummina short RNA-seq paired-end reads are first aligned on the consensus sequence of a TE family and then on the reference genome. The goal is to identify the number of fragments that are completely aligned inside the TEs ("Inside" fragments) and the number of fragments with a read aligned inside a TE and the other one mapping on the reference genome ("Outside" fragments). The rationale is that the sequencing of autonomously transcribed TEs should produce mostly paired-end reads mapping within the internal part of TE RNAs, hence mainly producing "Inside" fragments. Accordingly, by computing the Inside/Outside ratio, high ratio levels should result when TEs are autonomously transcribed.

# Technical steps:
1. Alignment of the reads (```--FASTQ_1 and --FASTQ_2 arguments```) on the TE consensus sequence (```--CONSENSUS argument```) with BWA
2. Not primary and supplementary alignments are discarded using Samtools 
3. The remaining fragments (referred to as a pair of reads) are filtered out if at least one read of the pair has more than a certain percentage of nucleotides that are not perfectly matched on the TE consensus sequence (```--MISMATCHES argument```)
4. The fragments that are completely aligned inside the TE consensus sequence are labeled as “Inside” fragments while the fragments with exclusively one read aligned inside the TE consensus sequence are labeled as “Outside” fragments 
5. Alignment of these sets of fragments on the reference genome with BWA
6. Application of the same previously described filters on the resulting BAM file
7. Facultative intersection (```--BED argument```) of the mapping coordinates of the reads with a BED file containing the TEs annotated on the reference genome (e.g. RepeatMasker track retrieved from UCSC Table Browser); fragments with at least one read aligned on an annotated TE are kept into account for the final calculation
8. The ratio between the number of Inside and Outside fragments is used as an indicator for the TE autonomous transcription level

# Prerequisites
* R packages: data.table, seqinr, here, stringr, ggplot2, tidyr, dplyr, magrittr, argparser, Rsamtools and GenomicAlignments
* Software: Samtools, Bedtools, BWA

# How to run Rindalyzer
The tool can be used from the command-line by running: 

```
cd $Rindalyzer
./Rindalyzer.R [arguments]
```

The arguments required are:
```
-1 / --FASTQ_1, Path to FASTQ r1 file (.fastq, .fq or .gz)
-2 / --FASTQ_2, Path to FASTQ r2 file (.fastq, .fq or .gz)
-c / --CONSENSUS, Path to FASTA file of consensus sequence (.fasta, .fa or .gz)
-m / --MISMATCHES, Maximum percentage of allowed mismatches (See Technical step 3)
-g / --GENOME, Path to FASTA file of the reference genome (.fasta, .fa or .gz)
-o / --OUTPUT, Path to directory for the output
-n / --NAME, Name of the analysis
```
The optional arguments are:
```
-t / --THREADS, Number of threads [default: 1]
-b / --BED, Path to BED file of genomic regions [default: none] (See Technical step 7)
-k / --KEEP_TEMPORARY, TRUE for maintaining temporary files [default: FALSE]
-s / --SAMTOOLS, Path to directory containing the Samtools software [default: none]; use if the path is not in the .bashrc
-B / --BEDTOOLS, Path to directory containing the Bedtools software [default: none]; use if the path is not in the .bashrc
-w / --BWA, Path to directory containing the BWA software [default: none]; use if the path is not in the .bashrc
```
 
Example usage:
```
./Rindalyzer.R --FASTQ_1 /path/to/fastq1.fastq --FASTQ_2 /path/to/fastq2.fastq --CONSENSUS /path/to/consensus.fasta --BED /path/to/file.bed --MISMATCHES 10 --GENOME /path/to/genome.fasta --OUTPUT /path/to/output/ --NAME test [ … ]
```

# Output files
1. ```rindReport.html``` is an HTML report for the analysis (the "GenomePostFilter" step refers to the final step of the analysis)
2. ```rindResults.txt``` is a tabular form of the fragment counts
3. ```Rindalyzer.log``` is a log file

# Development and Help

The Rindalyzer tool is developed by Mauro Esposito, PhD student in the Computational Genomics lab (SISSA/ISAS - Trieste - Italy) of prof. Remo Sanges. Please feel free to report bugs or suggestions.

# Citation

If you find Rindalyzer useful for your research, please cite: Esposito, M. et al. Transposons Acting as Competitive Endogenous RNAs: In-Silico Evidence from Datasets Characterised by L1 Overexpression. Biomedicines 2022, 10, 3279. https://doi.org/10.3390/biomedicines10123279
