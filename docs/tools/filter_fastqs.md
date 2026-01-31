# spatula filter-fastqs

## Summary 

`spatula filter-fastqs` filters a pair of FASTQ files based on whether the reads match specific IUPAC patterns. It can either keep or remove reads that satisfy the matching criteria, allowing for a specified number of mismatches.

## Required Options

* `--fq1 STR` : Input FASTQ file for read 1.
* `--fq2 STR` : Input FASTQ file for read 2.
* `--out1 STR` : Output FASTQ file for read 1.
* `--out2 STR` : Output FASTQ file for read 2.

## Additional Options

* `--pat1 STR` : IUPAC pattern to match against Read 1. Can be specified multiple times. At least one pattern must match if provided.
* `--pat2 STR` : IUPAC pattern to match against Read 2. Can be specified multiple times. At least one pattern must match if provided.
* `--min-mismatch INT` : Minimum number of mismatches allowed for a pattern to be considered a match. (Default: 0)
* `--remove` : If enabled, the tool removes the matching sequences from the output instead of keeping them. (Default: Keep matching sequences)

## Expected Output

The tool produces two output FASTQ files (specified by `--out1` and `--out2`) containing the filtered reads. 
* If `--remove` is OFF (default), only read pairs that match the specified patterns (according to `--pat1` and/or `--pat2`) are written to the output.
* If `--remove` is ON, read pairs that match the patterns are excluded, and only non-matching read pairs are written to the output.
* If no patterns are provided, all reads are preserved (unless `--remove` is used, in which case all reads are removed).
* The output files will be compressed if the filename ends with `.gz`.

## Full Usage 

The full usage of `spatula filter-fastqs` can be viewed with the `--help` option:

```
$ ./spatula filter-fastqs --help         
[./spatula filter-fastqs] -- Filter FASTQs based on a pattern

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --fq1          [STR: ]             : FASTQ file for read 1
   --fq2          [STR: ]             : FASTQ file for read 2

== Filtering options ==
   --pat1         [V_STR: ]           : IUPAC pattern to match for Read 1 (1+ match required)
   --pat2         [V_STR: ]           : IUPAC pattern to match for Read 2 (1+ match required)
   --min-mismatch [INT: 0]            : Minimum number of mismatches allowed per pattern
   --remove       [FLG: OFF]          : Remove the matching sequence from the FASTQ file instead of keeping it

== Output Options ==
   --out1         [STR: ]             : Output FASTQ file for read1
   --out2         [STR: ]             : Output FASTQ file for read2


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```