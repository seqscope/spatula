# spatula split-molecule-counts

!!! warning "Deprecated"
    `spatula split-molecule-counts` is **deprecated** and retained only for backward compatibility. Its functionality has been split into two separate commands:

    * [assign-feature2bin](assign_feature2bin.md) — assigns genes to bins from a feature count TSV and writes the bin assignment JSON.
    * [split-mol2bin](split_mol2bin.md) — splits a molecule TSV into per-bin files using that JSON.

    New workflows should use those two commands instead.

## Summary

`spatula split-molecule-counts` performs two tasks in a single command: (1) it reads a feature count TSV (`--feature-tsv`) and assigns each gene into ordered bins holding a roughly equal number of molecules, writing the assignment as a `_bin_counts.json` file; and (2) it uses that assignment to split the molecule TSV (`--mol-tsv`) into per-bin molecule and feature files, along with an index file.

Because the two steps are now available independently as [assign-feature2bin](assign_feature2bin.md) and [split-mol2bin](split_mol2bin.md), this combined command is deprecated.

## Required Options

* `--mol-tsv STR` : Input TSV file containing individual molecules.
* `--feature-tsv STR` : Input TSV file containing gene-level total counts (gene name in column 1, count in column 2).
* `--out-prefix STR` : Output prefix for the output TSV/JSON files.

## Additional Options

### Key Parameters
* `--bin-count INT` : Number of bins to split the features into (roughly equal number of molecules per bin). (Default: 50)
* `--skip-original` : Skip writing the original (unsplit) `_all_` files.
* `--compact-bin` : Write only the minimal columns (`X`, `Y`, gene, count) to the per-bin molecule files.

### Expected Columns in Input and Output
* `--colname-feature STR` : Column name for the gene name. (Default: `gene`)
* `--colname-count STR` : Column name for the gene count. (Default: `count`)
* `--colname-x STR` : Column name for the X coordinate. (Default: `X`)
* `--colname-y STR` : Column name for the Y coordinate. (Default: `Y`)
* `--col-rename V_STR` : Columns to rename in the output. Format: `old_name1:new_name1 old_name2:new_name2 ...`. May be specified multiple times.

### Auxiliary Input/Output Parameters
* `--in-mol-tsv-delim STR` : Delimiter for the input molecule TSV file. (Default: `\t`)
* `--in-feature-tsv-delim STR` : Delimiter for the input feature TSV file. (Default: `\t`)
* `--out-mol-tsv-delim STR` : Delimiter for the output molecule TSV files. (Default: `\t`)
* `--out-feature-tsv-delim STR` : Delimiter for the output feature TSV files. (Default: `\t`)
* `--out-mol-suffix STR` : Suffix for the output molecule TSV files. (Default: `molecules.tsv.gz`)
* `--out-json-suffix STR` : Suffix for the output JSON file containing bin counts and other metadata. (Default: `_bin_counts.json`)
* `--out-feature-suffix STR` : Suffix for the output feature TSV files. (Default: `features.tsv.gz`)
* `--strip-comment-char STR` : Character that marks comment lines to be skipped in the input files. (Default: `#`)

## Expected Output

Given `N` bins and output prefix `[out-prefix]`:

* `[out-prefix]_bin[i]_molecules.tsv.gz` and `[out-prefix]_bin[i]_features.tsv.gz` : Per-bin molecule and feature files for bin `i` (1-based).
* `[out-prefix]_all_molecules.tsv.gz` and `[out-prefix]_all_features.tsv.gz` : Unsplit files. Not written when `--skip-original` is set.
* `[out-prefix]_index.tsv` : Index file with one row per bin.
* `[out-prefix]_bin_counts.json` : JSON file with the gene → bin assignment (`{"gene", "count", "bin"}` per gene; `bin` is 1-based, `0` = unassigned).

## Full Usage

The full usage of `spatula split-molecule-counts` can be viewed with the `--help` option:

```
$ ./spatula split-molecule-counts --help
[./spatula split-molecule-counts] -- (Deprecated: use assign-feature2bin + split-mol2bin) Split molecule counts into pixel-level bins based on gene names

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --mol-tsv               [STR: ]             : TSV file containing individual molecules
   --feature-tsv           [STR: ]             : TSV file containing gene-level total counts (not used in the current version but can be used for filtering lowly expressed genes in the future)
   --out-prefix            [STR: ]             : Output Prefix for the joined TSV files

== Key Parameters ==
   --bin-count             [INT: 50]           : When --equal-bins is used, determine the number of bins to split the data into the same bin
   --skip-original         [FLG: OFF]          : Whether to skip writing the original file
   --compact-bin           [FLG: OFF]          : Compact the bins to store minimal information

== Expected columns in input and output ==
   --colname-feature       [STR: gene]         : Column name for gene name
   --colname-count         [STR: count]        : Column name for gene count
   --colname-x             [STR: X]            : Column name for X coordinate
   --colname-y             [STR: Y]            : Column name for Y coordinate
   --col-rename            [V_STR: ]           : Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...

== Auxilary Input/Output Parameters ==
   --in-mol-tsv-delim      [STR: 	]            : Delimiter for the input molecule TSV file
   --in-feature-tsv-delim  [STR: 	]            : Delimiter for the input feature TSV file
   --out-mol-tsv-delim     [STR: 	]            : Delimiter for the output molecule TSV file
   --out-feature-tsv-delim [STR: 	]            : Delimiter for the output feature TSV file
   --out-mol-suffix        [STR: molecules.tsv.gz] : Suffix for the output molecule TSV file
   --out-json-suffix       [STR: _bin_counts.json] : Suffix for the output JSON file containing bin counts and other metadata
   --out-feature-suffix    [STR: features.tsv.gz] : Suffix for the output feature TSV file
   --strip-comment-char    [STR: #]            : Character to strip from the beginning of lines in the input files (if any)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
