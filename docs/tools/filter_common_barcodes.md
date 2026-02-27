# spatula filter-common-barcodes

## Summary 

`spatula filter-common-barcodes` filters a pair of FASTQ files by retaining only the reads that contain a combination of Spatial Barcode (SBCD) and Unique Molecular Identifier (UMI) that appears at least a certain number of times (specified by `--min-freq`) across the entire dataset. This is useful for filtering out singleton or rare barcode-UMI combinations that might be errors.

## Required Options

* `--fq1 STR` : Input FASTQ file for read 1.
* `--fq2 STR` : Input FASTQ file for read 2.
* `--out STR` : Output prefix for the filtered FASTQ files.

## Additional Options

* `--start-sbcd INT` : 1-based start position of the spatial barcode (SBCD). (Default: 1)
* `--len-sbcd INT` : Length of the spatial barcode to read. (Default: 27)
* `--fq2-sbcd` : If set, the spatial barcode is extracted from Read 2. (Default: Read 1)
* `--start-umi INT` : 1-based start position of the UMI. (Default: 31)
* `--len-umi INT` : Length of the UMI to read. (Default: 9)
* `--fq2-umi` : If set, the UMI is extracted from Read 2. (Default: Read 1)
* `--min-freq INT` : Minimum frequency of SBCD-UMI pairs required to be considered valid and retained in the output. (Default: 2)
* `--nbatch INT` : Number of temporary batch files to use for processing. (Default: 512)

## Expected Output

The tool generates two filtered FASTQ files based on the specified output prefix:
* `[out_prefix].filt.R1.fastq.gz` : Filtered FASTQ file for Read 1.
* `[out_prefix].filt.R2.fastq.gz` : Filtered FASTQ file for Read 2.

These files contain only the read pairs where the combined SBCD and UMI sequence occurred at least `--min-freq` times in the input data.

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