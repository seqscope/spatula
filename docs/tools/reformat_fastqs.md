# spatula reformat-fastqs

## Summary 

`spatula reformat-fastqs` takes paired-end FASTQ files and reformats it into a format that can be used by `STARsolo` by filtering out irrelevant reads while arranging the barcodes and UMI sequences. Specifically, the software achieves the following:
* Filter out reads that does not have matching barcode sequences in the spatial barcode dictionary.
* Reformat the Read 1 sequences to be compatible with `STARsolo`, by copying the UMI sequences next to the spatial barcode sequences in Read 1.
* Trim Read 2 sequences to a specified length, which typically helps increase the mapping rates.
* If needed, trim the beginning of Read 1 sequences to match the spatial barcode sequences.

Here is a summary of a typical use case:

* Input: Takes (a) a pair of FASTQ files and (b) tsv files containing spatial barcodes to be matched.
* Output: Produces a pair of (uncompressed) FASTQ files after filtering reads and rearranging the part of sequences. 

A typical example is as follows:

```sh
spatula reformat-fastqs --fq1 /path/to/input.R1.fastq.gz \
                        --fq2 /path/to/input.R2.fastq.gz \
                        --match-tsv /path/to/match_sbcds_prefix1.match.sorted.uniq.tsv.gz \
                        --match-tsv /path/to/match_sbcds_prefix2.match.sorted.uniq.tsv.gz \
                        --match-tsv /path/to/match_sbcds_prefix3.match.sorted.uniq.tsv.gz \
                        --skip-sbcd 1 \
                        --len-sbcd 31 \
                        --len-umi 9 \
                        --len-r2 101 \
                        --out1 /path/to/output.R1.fastq \
                        --out2 /path/to/output.R2.fastq
```

See below for a more detailed usage description.

## Required options
* `--fq1` : The path to the FASTQ file of Read 1, which typically contains spatial barcode sequences.
* `--fq2` : The path to the FASTQ file of Read 2, which typically contains UMI sequences and cDNA sequences.
* `--match-tsv` : (Multiples allowed) The output `match.sorted.uniq.tsv.gz` file from `spatula match-sbcds` that contains the spatial barcode sequences to be matched. 
* `--out1` : The path to the output FASTQ files of Read 1.
* `--out2` : The path to the output FASTQ files of Read 2. See [Expected Output](#expected-output) for more details.

## Additional Options

* `--len-match` : The length of the spatial barcode to consider for a match. The default is 27, and the maximum possible value is 27.
* `--len-sbcd` : The length of the spatial barcode sequence to be copied in Read 1. The default is 30, and must be equal or smaller than actual spatial barcode length.
* `--len-umi` : The length of the UMI sequence (randomer) to be copied from Read 2 (beginning) to Read 1 (after spatial barcode). The default is 9.
* `--len-r2` : The length of Read 2 sequences to be trimmed. The default is 101.
* `--skip-sbcd` : The number of bases to be skipped in the beginning of the read. This is useful when insufficient bases are sequenced in the 1st-seq spatial barcode, which is reverse complemented. The default is 1.

## Expected Output

* `[out1]` is output FASTQ file for Read 1. It should contain spatial barcodes, followed by the UMI (randomer) sequences.
* `[out2]` is output FASTQ file for Read 2. It should contain the cDNA sequences, which includes the randomer sequence, trimmed to the specified length.

## Full Usage 

The full usage of `spatula reformat-fastqs` can be viewed with the `--help` option:

```
$ ./spatula reformat-fastqs --help
[./spatula reformat-fastqs] -- Reformat FASTQs to be ready for STARsolo alignment

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --fq1       [STR: ]             : FASTQ file for read 1
   --fq2       [STR: ]             : FASTQ file for read 2

== Filtering options ==
   --match-tsv [V_STR: ]           : .matched.tsv.gz file(s) to filter FASTQ file based on spatial barcodes
   --len-match [INT: 27]           : Length of perfect match required with 2nd-seq FASTQ (27 or less)

== Output Options ==
   --out1      [STR: ]             : Output FASTQ file for read1
   --out2      [STR: ]             : Output FASTQ file for read2
   --skip-sbcd [INT: 0]            : Skip first bases of spatial barcode (Read 1)
   --len-sbcd  [INT: 30]           : Length of spatial barcode (Read 1)
   --len-umi   [INT: 9]            : Length of UMI or randomer (Read 2)
   --len-r2    [INT: 101]          : Length of Read 2 to trim (including randomers)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```