# spatula split-mol2bin

## Summary

`spatula split-mol2bin` splits a molecule-level TSV file into multiple per-bin files based on a precomputed gene → bin assignment. The assignment is provided as a JSON file (the `--bin-json` file, produced by [assign-feature2bin](assign_feature2bin.md) and equivalent to the old `_bin_counts.json`). Each molecule is routed to the bin of its gene, and per-bin molecule TSV files, per-bin feature (gene count) TSV files, and an index file are written.

This is the second of the two commands that replace the (deprecated) [split-molecule-counts](split_molecule_counts.md) command. It performs only the splitting step; the gene → bin assignment is computed separately by [assign-feature2bin](assign_feature2bin.md).

## Required Options

* `--mol-tsv STR` : Input TSV file containing individual molecules (with per-molecule coordinates, gene, and count columns).
* `--bin-json STR` : JSON file containing the gene → bin assignment produced by [assign-feature2bin](assign_feature2bin.md) (equivalent to the old `_bin_counts.json`).
* `--out-prefix STR` : Output prefix for the per-bin TSV files.

## Additional Options

### Key Parameters
* `--skip-original` : Skip writing the original (unsplit) `_all_` molecule and feature files, only writing the per-bin files.
* `--compact-bin` : Write only the minimal columns (`X`, `Y`, gene, count) to the per-bin molecule files instead of copying all input columns.

### Expected Columns in Input and Output
* `--colname-feature STR` : Column name for the gene name. (Default: `gene`)
* `--colname-count STR` : Column name for the gene count. (Default: `count`)
* `--colname-x STR` : Column name for the X coordinate. (Default: `X`)
* `--colname-y STR` : Column name for the Y coordinate. (Default: `Y`)
* `--col-rename V_STR` : Columns to rename in the output. Format: `old_name1:new_name1 old_name2:new_name2 ...`. May be specified multiple times.

### Auxiliary Input/Output Parameters
* `--in-mol-tsv-delim STR` : Delimiter for the input molecule TSV file. (Default: `\t`)
* `--out-mol-tsv-delim STR` : Delimiter for the output molecule TSV files. (Default: `\t`)
* `--out-feature-tsv-delim STR` : Delimiter for the output feature TSV files. (Default: `\t`)
* `--out-mol-suffix STR` : Suffix for the output molecule TSV files. (Default: `molecules.tsv.gz`)
* `--out-feature-suffix STR` : Suffix for the output feature TSV files. (Default: `features.tsv.gz`)
* `--strip-comment-char STR` : Character that marks comment lines to be skipped in the input files. (Default: `#`)

## Expected Output

Given `N` bins (as determined by the `--bin-json` assignment) and output prefix `[out-prefix]`:

* `[out-prefix]_bin[i]_[out-mol-suffix]` : Per-bin molecule TSV file for bin `i` (1-based).
* `[out-prefix]_bin[i]_[out-feature-suffix]` : Per-bin feature (gene, count) TSV file for bin `i`.
* `[out-prefix]_all_[out-mol-suffix]` : All molecules (unsplit). Not written when `--skip-original` is set.
* `[out-prefix]_all_[out-feature-suffix]` : All features (unsplit). Not written when `--skip-original` is set.
* `[out-prefix]_index.tsv` : Index file with one row per bin (and an `all` row unless `--skip-original`), listing `bin_id`, `molecule_count`, `features_count`, `molecules_path`, and `features_path`.

## Full Usage

The full usage of `spatula split-mol2bin` can be viewed with the `--help` option:

```
$ ./spatula split-mol2bin --help
[./spatula split-mol2bin] -- Split a molecule TSV into per-bin files using a gene->bin assignment JSON (from assign-feature2bin)

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --mol-tsv               [STR: ]             : TSV file containing individual molecules
   --bin-json              [STR: ]             : JSON file containing the gene->bin assignment (produced by assign-feature2bin, equivalent to the old _bin_counts.json)
   --out-prefix            [STR: ]             : Output Prefix for the per-bin TSV files

== Key Parameters ==
   --skip-original         [FLG: OFF]          : Whether to skip writing the original (unsplit) file
   --compact-bin           [FLG: OFF]          : Compact the bins to store minimal information

== Expected columns in input and output ==
   --colname-feature       [STR: gene]         : Column name for gene name
   --colname-count         [STR: count]        : Column name for gene count
   --colname-x             [STR: X]            : Column name for X coordinate
   --colname-y             [STR: Y]            : Column name for Y coordinate
   --col-rename            [V_STR: ]           : Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...

== Auxilary Input/Output Parameters ==
   --in-mol-tsv-delim      [STR: 	]            : Delimiter for the input molecule TSV file
   --out-mol-tsv-delim     [STR: 	]            : Delimiter for the output molecule TSV file
   --out-feature-tsv-delim [STR: 	]            : Delimiter for the output feature TSV file
   --out-mol-suffix        [STR: molecules.tsv.gz] : Suffix for the output molecule TSV file
   --out-feature-suffix    [STR: features.tsv.gz] : Suffix for the output feature TSV file
   --strip-comment-char    [STR: #]            : Character to strip from the beginning of lines in the input files (if any)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
