# crispr-parsr
Software to parse and analyse deletions (and insertions) from a CRISPR resequencing experiment

## Synopsis

crispr-parsr.pl [--wt ref_seq.fa] [--samples merge.txt] --input INPUT_DIR --output OUTPUT_DIR

## Description

Briefly, the pipeline does the following:

1. Index the wild-type sequence for Bowtie2
2. Merge FASTQ files according to the information in the samples file
3. Trim and QC files with Trim Galore! (which uses cutadapt and FastQC internally)
4. Align the reads with Bowtie2
5. Parses the alignment to extract the deletions (and insertions) found in each sample

## Options

**--input INPUT_DIR**

The directory where all the FASTQ files are. By default, crispr-parsr will look for a FASTA file
(for the wild-type sequence) and a TXT file (for the samples definition) in this directory as well.

**--output OUTPUT_DIR**

All files with intermediary data and final plots will be created in this directory. There will be a
README.txt file explaining the content of each file.

**--wt ref_seq.fa**

This is a simple FASTA file with one sequence only (typically a short one) corresponding to the expected
wild-type sequence, before editing has ocurred.

If the filename is not provided, crispr-parsr will look for a FASTA file within the INPUT_DIR. If either
none or more than one file with the .fa extension exists in the INPUT_DIR, this will fail.

You can either provide the full path to the file or simply the name of the file if it is located
in the INPUT_DIR.

Example:
```
>sample_amplicon1_ref_123456
AACAGTGGGTCTTCGGGAACCAAGCAGCAGCCATGGGAGGTCTGGCTGTGTTCAGGCTCTGCTCGTGTAGATTCACAGCGCGCTCTGAACCCCCGCTG
AGCTACCGATGGAAGAGGAGGAGGTCCTACAGTCGGAGATTCACAGCGCGCTCTGAACCACTTTCAGGAGACTCGACTATTATGACTTATACGCGATA
```

**--samples merge.txt**

This file contains the relation of FASTQ files for each sample. As the pipeline has been designed for
paired end reads, you need to specify the list of R1 and R2 fastq files in consecutive lines. You can
have more than one set of PE reads for each sample. For this, you need to separate your filenames with
a [tab] (you can use Excel and save the file in tab-separated format). Please use the same order for
R1 and R2 files within a sample.

You can also specify a label for each sample (recommended).

If the filename is not provided, crispr-parsr will look for a TXT file within the INPUT_DIR. If either
none or more than one file with the .txt extension exists in the INPUT_DIR, this will fail.

You can either provide the full path to the file or simply the name of the file if it is located
in the INPUT_DIR.

Example (with labels):
```
Mock:
12345A01_S12_L001_R1_001.fastq
12345A01_S12_L001_R2_001.fastq
Run1:
12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
Run2:
12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
```

Example (without labels):
```
12345A01_S12_L001_R1_001.fastq
12345A01_S12_L001_R2_001.fastq
12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
```
