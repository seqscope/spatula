# spatula subset-mex

## Summary

`spatula subset-mex` subsets a 10x Market Exchange (MEX) formatted digital gene expression (DGE) matrix to a specific set of features and/or barcodes. Features and barcodes can be filtered by inclusion/exclusion lists and by minimum total counts. The output is written as a new MEX triplet (barcodes, features, matrix) in the specified output directory.

See [mex2sptsv](mex2sptsv.md) and [sptsv2mex](sptsv2mex.md) for converting between MEX and sparse TSV formats.

## Required Options

* `--in-dir STR` : Input directory containing the MEX files.
* `--out-dir STR` : Output directory for the subset MEX files.

## Additional Options

### Input Files
* `--in-bcd STR` : Input barcode file name. (Default: `barcodes.tsv.gz`)
* `--in-ftr STR` : Input feature file name. (Default: `features.tsv.gz`)
* `--in-mtx STR` : Input matrix file name. (Default: `matrix.mtx.gz`)
* `--icol-mtx INT` : 1-based column index in the matrix file to use as the count. (Default: 3)

### Output Files
* `--out-bcd STR` : Output barcode file name. (Default: `barcodes.tsv.gz`)
* `--out-ftr STR` : Output feature file name. (Default: `features.tsv.gz`)
* `--out-mtx STR` : Output matrix file name. (Default: `matrix.mtx.gz`)

### Filtering
* `--include-feature-list STR` / `--exclude-feature-list STR` : File listing features (genes) to include/exclude.
* `--include-barcode-list STR` / `--exclude-barcode-list STR` : File listing barcodes to include/exclude.
* `--min-feature-count INT` : Minimum total count for a feature to be included. Requires reading the input matrix twice. (Default: 0)
* `--min-barcode-count INT` : Minimum total count for a barcode to be included. Requires reading the input matrix twice. (Default: 0)

## Expected Output

* `[out-dir]/[out-bcd]` : Subset barcode file.
* `[out-dir]/[out-ftr]` : Subset feature file.
* `[out-dir]/[out-mtx]` : Subset matrix file.

## Full Usage

The full usage of `spatula subset-mex` can be viewed with the `--help` option:

```
$ ./spatula subset-mex --help
[./spatula subset-mex] -- Subset 10x Market Exchange (MEX) formatted DGE to specific features and barcodes

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input Options ==
   --in-dir               [STR: ]             : Input directory containing the MEX files
   --in-bcd               [STR: barcodes.tsv.gz] : Input Barcode file name
   --in-ftr               [STR: features.tsv.gz] : Input Feature file name
   --in-mtx               [STR: matrix.mtx.gz] : Input Matrix file name
   --icol-mtx             [INT: 3]            : 1-based column index in the matrix file to use as the count

== Key Output Options ==
   --out-dir              [STR: ]             : Output directory containing the MEX files
   --out-bcd              [STR: barcodes.tsv.gz] : Output Barcode file name
   --out-ftr              [STR: features.tsv.gz] : Output Feature file name
   --out-mtx              [STR: matrix.mtx.gz] : Output Matrix file name

== Input Filtering Options ==
   --include-feature-list [STR: ]             : A file containing a list of input genes to be included (feature name of IDs)
   --exclude-feature-list [STR: ]             : A file containing a list of input genes to be excluded (feature name of IDs)
   --include-barcode-list [STR: ]             : A file containing a list of input barcode IDs to be included (feature name of IDs)
   --exclude-barcode-list [STR: ]             : A file containing a list of input barcode IDs to be excluded (feature name of IDs)
   --min-feature-count    [INT: 0]            : Minimum feature count to include in the output. Requires reading the input matrix twice
   --min-barcode-count    [INT: 0]            : Minimum barcode count to include in the output. Requires reading the input matrix twice


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
