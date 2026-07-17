# spatula stratify-fastq-by-barcodes

## Summary

`spatula stratify-fastq-by-barcodes` stratifies (groups) paired-end FASTQ reads by their spatial barcode (SBCD) and UMI. It extracts the spatial barcode and UMI from configurable positions in Read 1 or Read 2, and organizes reads into batches keyed by the SBCD–UMI pair. Pairs can be filtered by minimum and maximum frequency so that only barcode/UMI combinations within a chosen abundance range are retained.

## Required Options

* `--fq1 STR` : FASTQ file for Read 1.
* `--fq2 STR` : FASTQ file for Read 2.
* `--out STR` : Output prefix.

## Additional Options

### Barcode and UMI Positions
* `--start-sbcd INT` : 1-based start position of the spatial barcode. (Default: 1)
* `--len-sbcd INT` : Length of the spatial barcode to read. (Default: 27)
* `--fq2-sbcd` : The spatial barcode is located in Read 2 (instead of Read 1). (Default: OFF)
* `--start-umi INT` : 1-based start position of the UMI. (Default: 31)
* `--len-umi INT` : Length of the UMI to read. (Default: 9)
* `--fq2-umi` : The UMI is located in Read 2. (Default: OFF)

### Output / Filtering
* `--nbatch INT` : Number of batches to process. (Default: 512)
* `--min-freq INT` : Minimum frequency of an SBCD–UMI pair for it to be considered valid. (Default: 1)
* `--max-freq INT` : Maximum frequency of an SBCD–UMI pair for it to be considered valid. (Default: 2147483647)

## Expected Output

* `[out]` (prefix) : Batched output files of reads stratified by SBCD–UMI pair.

## Full Usage

The full usage of `spatula stratify-fastq-by-barcodes` can be viewed with the `--help` option:

```
$ ./spatula stratify-fastq-by-barcodes --help
[./spatula stratify-fastq-by-barcodes] -- Stratify FASTQ files by spatial barcodes and UMIs

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --fq1        [STR: ]             : FASTQ file for read 1
   --fq2        [STR: ]             : FASTQ file for read 2

== Barcode and UMI options ==
   --start-sbcd [INT: 1]            : 1-based start position of spatial barcode
   --len-sbcd   [INT: 27]           : Length of spatial barcode to read
   --fq2-sbcd   [FLG: OFF]          : Spatial barcode is in Read 2
   --start-umi  [INT: 31]           : 1-based start position of UMI
   --len-umi    [INT: 9]            : Length of UMI to read
   --fq2-umi    [FLG: OFF]          : UMI is in Read 2

== Output Options ==
   --out        [STR: ]             : Output prefix
   --nbatch     [INT: 512]          : Number of batches to process
   --min-freq   [INT: 1]            : Minimum frequency of SBCD-UMI pairs to be considered as a valid pair
   --max-freq   [INT: 2147483647]   : Maximum frequency of SBCD-UMI pairs to be considered as a valid pair


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
