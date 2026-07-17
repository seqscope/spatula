# spatula make-sdge

## Summary

`spatula make-sdge` makes spatial DGE files from STARsolo output by mapping target barcodes (or a target FASTQ) to a reference spatial barcode dictionary using spaced k-mers. It is an alternative barcode-matching path used when constructing spatial gene expression from raw reads, supporting approximate matching within a configurable edit distance.

## Required Options

* `--ref STR` : Reference barcode dictionary.
* `--out STR` : Output file name.
* A target must be provided via either `--tgt` (target barcode file) or `--fq` (target FASTQ file).

## Additional Options

### Input
* `--tgt STR` : Target barcode file to map to the reference.
* `--fq STR` : Target FASTQ file to map to the reference (trims and ignores qualities).

### Spaced k-mer Parameters
* `--len INT` : Expected length of the barcode. (Default: 0)
* `--lunit INT` : Expected length of k-mer units (must be a multiple of 4). (Default: 0)
* `--nunit INT` : Expected number of k-mer units required to form a single key to declare a match. (Default: 0)
* `--chunk-size INT` : Size of chunks used to store spaced k-mers. (Default: 100000)
* `--max-dist INT` : Maximum edit distance to allow. (Default: -1)
* `--skip-no-match` : Skip reads that do not match.
* `--skip-xtra` : Skip writing extra fields even when available.

### Output
* `--verbose INT` : Verbosity threshold. `0` = fully verbose, `100` = no verbose messages. (Default: 100)

## Expected Output

* `[out]` : The spatial DGE output file mapping target barcodes/reads to the reference dictionary.

## Full Usage

The full usage of `spatula make-sdge` can be viewed with the `--help` option:

```
$ ./spatula make-sdge --help
[./spatula make-sdge] -- Make spatial DGE files from STARsolo output

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --ref           [STR: ]             : Reference barcode dictionary
   --tgt           [STR: ]             : Target barcode to map to the reference
   --fq            [STR: ]             : Target FASTQ file to map to the reference, trimming and ignoring qualities

== Parameters for spaced k-mers ==
   --len           [INT: 0]            : Expected length of barcode
   --lunit         [INT: 0]            : Expected length of kmer units (must be a multiple of 4)
   --nunit         [INT: 0]            : Expected number of kmer units to make a single key to declare a match
   --chunk-size    [INT: 100000]       : Size of chunks to store spaced kmers
   --max-dist      [INT: -1]           : Maximum edit distance to allow
   --skip-no-match [FLG: OFF]          : Skip reads that does not match
   --skip-xtra     [FLG: OFF]          : Skip writing extra fields even if it is available

== Output Options ==
   --out           [STR: ]             : Output file name
   --verbose       [INT: 100]          : Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
