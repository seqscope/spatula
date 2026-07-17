# spatula search-tag

## Summary

`spatula search-tag` searches for tag sequences in paired FASTQ files against a tag dictionary, using a whitelist of spatial barcodes from the `smatch` step. It is similar to [match-tag](match_tag.md) but is oriented toward searching/scanning for tags and does not take a separate tag-position argument (only spatial barcode and optional UMI positions are specified).

## Required Options

* `--fq1 STR` : FASTQ file Read 1 containing the 2nd-seq spatial barcode.
* `--fq2 STR` : FASTQ file Read 2 containing the tag (and optionally UMI).
* `--tag STR` : Dictionary file containing tag sequences to search for.
* `--smatch STR` : Output from the `smatch` step containing spatial barcodes used as the whitelist.
* `--out STR` : Output prefix to store the tags.

## Additional Options

### Positions in Reads
* `--bcd-pos STR` : 1-based positions of the spatial barcode sequences in Read 1, in `[beg1]-[end1],[beg2]-[end2],...` format. (Default: `1-32`)
* `--umi-pos STR` : 1-based positions of the UMI sequences in Read 2, in the same format.

### Batch
* `--batch INT` : Size of a single batch to store. (Default: 10000000)

## Expected Output

* `[out]` (prefix) : Output files reporting the tags found for each spatial barcode.

## Full Usage

The full usage of `spatula search-tag` can be viewed with the `--help` option:

```
$ ./spatula search-tag --help
[./spatula search-tag] -- Search tags from FASTQ

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input files ==
   --fq1     [STR: ]             : FASTQ file read 1 containing 2nd-seq spatial barcode
   --fq2     [STR: ]             : FASTQ file read 2 containing 2nd-seq spatial barcode
   --tag     [STR: ]             : Dictionary file containing tag sequences to search for
   --smatch  [STR: ]             : Output from smatch step that contains spatial barcodes to be used for whitelist

== Positions in reads ==
   --bcd-pos [STR: 1-32]         : String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of Spatial barcode sequences (in Read 1)
   --umi-pos [STR: ]             : String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of UMI sequences (in Read 2)

== Batch options ==
   --batch   [INT: 10000000]     : Size of single batch to store stored files

== Output Options ==
   --out     [STR: ]             : Output prefix to store the tags


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
