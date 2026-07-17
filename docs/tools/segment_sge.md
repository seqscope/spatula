# spatula segment-sge

## Summary

`spatula segment-sge` segments a spatial gene expression (SGE) matrix into spatial units using custom, grid, or hexagonal masks (supplied as NumPy `.npy` files). Transcripts are aggregated into the mask-defined segments, producing a new SGE where each unit corresponds to a segment rather than an individual barcode.

## Required Options

* `--in-sge STR` : Input SGE directory.
* `--in-npy STR` : Input SGE `.npy` mask file(s) defining the segmentation.
* `--out-sge STR` : Prefix of the output SGE directory.

## Additional Options

### SGE Input
* `--bcd STR` : Barcode file name. (Default: `barcodes.tsv.gz`)
* `--ftr STR` : Feature file name. (Default: `features.tsv.gz`)
* `--mtx STR` : Matrix file name. (Default: `matrix.mtx.gz`)
* `--icol-mtx INT` : 1-based column index in the matrix file to use as the count. (Default: 1)
* `--icol-bcd-barcode INT` : 1-based column index of the barcode in the barcode file. (Default: 1)
* `--icol-bcd-x INT` / `--icol-bcd-y INT` : 1-based column index of the X/Y coordinate in the barcode file. (Defaults: 6, 7)
* `--icol-ftr-id INT` / `--icol-ftr-name INT` : 1-based column index of the feature ID/name in the feature file. (Defaults: 1, 2)

### Output
* `--units-per-px FLT` : Coordinate units per pixel. (Default: `1.00`)

## Expected Output

* `[out-sge]` : Output SGE directory in which each spatial unit corresponds to a mask segment.

## Full Usage

The full usage of `spatula segment-sge` can be viewed with the `--help` option:

```
$ ./spatula segment-sge --help
[./spatula segment-sge] -- Segment SGE based on custom, grid, or hexagonal masks

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --in-sge           [STR: ]             : Input SGE directory
   --in-npy           [STR: ]             : Input SGE npy files
   --out-sge          [STR: ]             : Prefix of output SGE directory

== Input Options for SGE input ==
   --bcd              [STR: barcodes.tsv.gz] : Barcode file name
   --ftr              [STR: features.tsv.gz] : Feature file name
   --mtx              [STR: matrix.mtx.gz] : Matrix file name
   --icol-mtx         [INT: 1]            : 1-based column index in the matrix file to use as the count
   --icol-bcd-barcode [INT: 1]            : 1-based column index of barcode in the barcode file
   --icol-bcd-x       [INT: 6]            : 1-based column index of x coordinate in the barcode file
   --icol-bcd-y       [INT: 7]            : 1-based column index of y coordinate in the barcode file
   --icol-ftr-id      [INT: 1]            : 1-based column index of feature ID in the barcode file
   --icol-ftr-name    [INT: 2]            : 1-based column index of feature name in the barcode file

== Key Output Options ==
   --units-per-px     [FLT: 1.00]         : Coordinate unit per pixel


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
