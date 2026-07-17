# spatula match-tag

## Summary

`spatula match-tag` matches tag sequences (e.g. antibody-derived tags, hashing tags, or other custom barcodes) from paired FASTQ files against a tag dictionary, using a whitelist of spatial barcodes from the `smatch` step. For each read pair it extracts the spatial barcode (from Read 1), optionally a UMI (from Read 2), and the tag sequence (from Read 2), and produces per-tag matched output that can be merged with [merge-matched-tags](merge_matched_tags.md).

See also [search-tag](search_tag.md) for searching tags without a UMI.

## Required Options

* `--fq1 STR` : FASTQ file Read 1 containing the 2nd-seq spatial barcode.
* `--fq2 STR` : FASTQ file Read 2 containing the tag (and optionally UMI).
* `--tag STR` : Dictionary file containing tag sequences.
* `--smatch STR` : Output from the `smatch` step containing spatial barcodes used as the whitelist.
* `--out STR` : Output prefix to store the matched tags.

## Additional Options

### Positions in Reads
* `--bcd-pos STR` : 1-based positions of the spatial barcode sequences in Read 1, in `[beg1]-[end1],[beg2]-[end2],...` format. (Default: `1-32`)
* `--umi-pos STR` : 1-based positions of the UMI sequences in Read 2, in the same format.
* `--tag-pos STR` : 1-based positions of the tag sequences in Read 2, in the same format. (Default: `1-15`)

### Batch
* `--batch INT` : Size of a single batch to store. (Default: 10000000)

## Expected Output

* `[out]` (prefix) : Per-tag matched output files associating spatial barcodes (and UMIs) with tags.

## Full Usage

The full usage of `spatula match-tag` can be viewed with the `--help` option:

```
$ ./spatula match-tag --help
[./spatula match-tag] -- Match tags from FASTQ

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input files ==
   --fq1     [STR: ]             : FASTQ file read 1 containing 2nd-seq spatial barcode
   --fq2     [STR: ]             : FASTQ file read 2 containing 2nd-seq spatial barcode
   --tag     [STR: ]             : Dictionary file containing tag sequences
   --smatch  [STR: ]             : Output from smatch step that contains spatial barcodes to be used for whitelist

== Positions in reads ==
   --bcd-pos [STR: 1-32]         : String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of Spatial barcode sequences (in Read 1)
   --umi-pos [STR: ]             : String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of UMI sequences (in Read 2)
   --tag-pos [STR: 1-15]         : String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of tag sequences (in Read 2)

== Batch options ==
   --batch   [INT: 10000000]     : Size of single batch to store stored files

== Output Options ==
   --out     [STR: ]             : Output prefix to store the tags


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
