# spatula sge2tsv

## Summary

`spatula sge2tsv` converts a spatial gene expression (SGE) matrix — in the barcode/feature/matrix format produced by `dge2sge`/`make-sdge` — into a plain, generic transcript TSV file (unsorted). Each nonzero matrix entry becomes a row with `X`, `Y`, gene, and count columns, with coordinates converted to micrometers. Features can be filtered by inclusion/exclusion lists, substrings, or regular expressions.

This is a lighter-weight sibling of [convert-sge](convert_sge.md); the output is not sorted by spatial coordinate and typically needs a subsequent `sort` step for downstream FICTURE analysis.

## Required Options

* `--sge STR` : Input SGE directory.
* `--out STR` : Output prefix.

## Additional Options

### Input Files
* `--bcd STR` : Barcode file name. (Default: `barcodes.tsv.gz`)
* `--ftr STR` : Feature file name. (Default: `features.tsv.gz`)
* `--mtx STR` : Matrix file name. (Default: `matrix.mtx.gz`)
* `--icol-mtx STR` : 1-based column index in the matrix file to use as the count.

### Feature Filtering
* `--include-feature-list STR` / `--exclude-feature-list STR` : File listing genes to include/exclude.
* `--include-feature-substr STR` / `--exclude-feature-substr STR` : Substring of feature names to include/exclude.
* `--include-feature-regex STR` / `--exclude-feature-regex STR` : Regex of feature names to include/exclude.

### Output
* `--units-per-um FLT` : Coordinate units per micrometer (conversion factor). (Default: `1.00`)
* `--precision-um INT` : Output precision below the decimal point. (Default: 3)
* `--print-feature-id` : Print feature ID in addition to the feature name.
* `--colname-feature-name STR` : Column name for the feature/gene name. (Default: `gene`)
* `--colname-feature-id STR` : Column name for the feature/gene ID. (Default: `MoleculeID`)
* `--colname-count STR` : Column name for the count. (Default: `Count`)
* `--colname-x STR` / `--colname-y STR` : Column names for X/Y. (Defaults: `X`, `Y`)
* `--suffix-tsv STR` : Suffix for the transcript output file. (Default: `.transcripts.tsv.gz`)
* `--suffix-ftr STR` : Suffix for the feature output file. (Default: `.features.clean.tsv.gz`)
* `--suffix-minmax STR` : Suffix for the minmax output file. (Default: `.mimmax.tsv`)

## Expected Output

* `[out][suffix-tsv]` : The generic transcript TSV file (unsorted).
* `[out][suffix-ftr]` : The cleaned feature list with total counts.
* `[out][suffix-minmax]` : The minmax coordinate file.

## Full Usage

The full usage of `spatula sge2tsv` can be viewed with the `--help` option:

```
$ ./spatula sge2tsv --help
[./spatula sge2tsv] -- Convert SGE (from sge2sdge) into plain TSV format (unsorted)

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input Options ==
   --sge                    [STR: ]             : SGE directory
   --bcd                    [STR: barcodes.tsv.gz] : Barcode file name
   --ftr                    [STR: features.tsv.gz] : Feature file name
   --mtx                    [STR: matrix.mtx.gz] : Matrix file name
   --icol-mtx               [STR:              : 1-based column index in the matrix file to use as the count

== Input Filtering Options ==
   --include-feature-list   [STR: ]             : A file containing a list of input genes to be included (feature name of IDs)
   --exclude-feature-list   [STR: ]             : A file containing a list of input genes to be excluded (feature name of IDs)
   --include-feature-substr [STR: ]             : A substring of feature/gene names to be included
   --exclude-feature-substr [STR: ]             : A substring of feature/gene names to be excluded
   --include-feature-regex  [STR: ]             : A regex pattern of feature/gene names to be included
   --exclude-feature-regex  [STR: ]             : A regex pattern of feature/gene names to be excluded

== Key Output Options ==
   --out                    [STR: ]             : Output prefix
   --units-per-um           [FLT: 1.00]         : Coordinate unit per um
   --precision-um           [INT: 3]            : Output precision below the decimal point

== Auxilary Output Options ==
   --print-feature-id       [FLG: OFF]          : Print feature ID in output file
   --colname-feature-name   [STR: gene]         : Column name for feature/gene name
   --colname-feature-id     [STR: MoleculeID]   : Column name for feature/gene ID
   --colname-count          [STR: Count]        : Column name for Count
   --colname-x              [STR: X]            : Column name for X
   --colname-y              [STR: Y]            : Column name for Y
   --suffix-tsv             [STR: .transcripts.tsv.gz] : Suffix for the transcript output file
   --suffix-ftr             [STR: .features.clean.tsv.gz] : Suffix for the feature output file
   --suffix-minmax          [STR: .mimmax.tsv]  : Suffix for the minmax output file


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
