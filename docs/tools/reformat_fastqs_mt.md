# spatula reformat-fastqs-mt

!!! warning "Experimental / not working"
    `spatula reformat-fastqs-mt` is an experimental multithreaded variant of [reformat-fastqs](reformat_fastqs.md) and is currently marked as **not working** in its own help text. Use [reformat-fastqs](reformat_fastqs.md) for production workflows.

## Summary

`spatula reformat-fastqs-mt` is intended to be a multithreaded version of [reformat-fastqs](reformat_fastqs.md), reformatting paired-end FASTQ files for `STARsolo` alignment: filtering reads whose spatial barcodes are not present in the match TSV, copying UMI sequences next to the spatial barcode in Read 1, and trimming Read 2. Unlike [reformat-fastqs](reformat_fastqs.md), the input and output FASTQ files are expected to be plain (uncompressed).

## Required Options

* `--fq1 STR` : FASTQ file for Read 1 (plain file).
* `--fq2 STR` : FASTQ file for Read 2 (plain file).
* `--match-tsv V_STR` : `.matched.tsv.gz` file(s) used to filter the FASTQ based on spatial barcodes. May be specified multiple times.
* `--out1 STR` : Output FASTQ file for Read 1 (plain file).
* `--out2 STR` : Output FASTQ file for Read 2 (plain file).

## Additional Options

* `--len-match INT` : Length of perfect match required with the 2nd-seq FASTQ (27 or less). (Default: 27)
* `--len-sbcd INT` : Length of the spatial barcode (Read 1). (Default: 30)
* `--len-umi INT` : Length of the UMI or randomer (Read 2). (Default: 9)
* `--len-r2 INT` : Length of Read 2 to trim (including randomers). (Default: 101)

## Expected Output

* `[out1]` / `[out2]` : Reformatted plain FASTQ files for Read 1 and Read 2.

## Full Usage

The full usage of `spatula reformat-fastqs-mt` can be viewed with the `--help` option:

```
$ ./spatula reformat-fastqs-mt --help
[./spatula reformat-fastqs-mt] -- Reformat FASTQs to be ready for STARsolo alignment (multithreaded/not working)

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --fq1       [STR: ]             : FASTQ file for read 1 (plain file)
   --fq2       [STR: ]             : FASTQ file for read 2 (plain file)

== Filtering options ==
   --match-tsv [V_STR: ]           : .matched.tsv.gz file(s) to filter FASTQ file based on spatial barcodes
   --len-match [INT: 27]           : Length of perfect match required with 2nd-seq FASTQ (27 or less)

== Output Options ==
   --out1      [STR: ]             : Output FASTQ file for read1 (plain file)
   --out2      [STR: ]             : Output FASTQ file for read2 (plain file)
   --len-sbcd  [INT: 30]           : Length of spatial barcode (Read 1)
   --len-umi   [INT: 9]            : Length of UMI or randomer (Read 2)
   --len-r2    [INT: 101]          : Length of Read 2 to trim (including randomers)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
