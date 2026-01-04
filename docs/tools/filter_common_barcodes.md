# spatula filter-common-barcodes

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

## Full Usage 

The full usage of `spatula filter-common-barcodes` can be viewed with the `--help` option:

```
./spatula filter-common-barcodes --help
[./spatula filter-common-barcodes] -- Filter FASTQ files based on the frequency of barcodes

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
   --min-freq   [INT: 2]            : Minimum frequency of SBCD-UMI pairs to be considered as a valid pair


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```