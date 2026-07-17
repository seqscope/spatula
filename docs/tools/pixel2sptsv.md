# spatula pixel2sptsv

## Summary

`spatula pixel2sptsv` converts a pixel-level (transcript-level) TSV into the sparse TSV (SpTSV) format used by FICTURE2, aggregating transcripts by a previously assigned cell ID. Each output unit corresponds to a cell, and its features are the summed transcript counts. The tool supports extensive filtering (by cell ID lists, feature lists, regex/substring, and minimum counts) and writes a metadata JSON, a cell metadata table, a feature count table, and the sparse TSV itself.

See [pixels2cells](pixels2cells.md) or [tsv-add-cell-id](tsv_add_cell_id.md) for producing the input cell IDs, and [sptsv2mex](sptsv2mex.md) for converting the resulting SpTSV to MEX.

## Required Options

* `--pixel STR` : Input pixel-level transcript data.
* `--out STR` : Output prefix.

## Additional Options

### Input Columns
* `--in-col-id STR` : Column name for cell ID. (Default: `cell_id`)
* `--in-col-ftr STR` : Column name for feature name. (Default: `gene`)
* `--in-col-cnt STR` : Column name for count. (Default: `count`)
* `--in-col-x STR` / `--in-col-y STR` : Column names for X/Y coordinates. (Defaults: `X`, `Y`)
* `--no-header` : Indicates the input file has no header line. When set, columns are referenced by 1-based index instead.
* `--idx-col-x INT` / `--idx-col-y INT` / `--idx-col-ftr INT` / `--idx-col-cnt INT` / `--idx-col-id INT` : 1-based column indices used only with `--no-header`. (Defaults: 1, 2, 3, 4, 5)

### Input Filtering
* `--ignore-ids STR` : Comma-separated IDs treated as null and ignored. (Default: `UNASSIGNED,-1,0,NA`)
* `--include-id-list STR` / `--exclude-id-list STR` : File listing input cell IDs to include/exclude.
* `--include-feature-list STR` / `--exclude-feature-list STR` : File listing genes to include/exclude.
* `--include-feature-regex STR` / `--exclude-feature-regex STR` : Regex of feature names to include/exclude.
* `--include-feature-substr STR` / `--exclude-feature-substr STR` : Substring of feature names to include/exclude.
* `--min-cell-count INT` : Minimum total count for a cell to be included. (Default: 1)
* `--min-feature-count INT` : Minimum total count for a feature to be included. (Default: 1)

### Output
* `--out-suffix-json STR` : Suffix for the metadata JSON file. (Default: `.json`)
* `--out-suffix-cell-meta STR` : Suffix for the cell metadata output. (Default: `.cell.metadata.tsv`)
* `--out-suffix-tsv STR` : Suffix for the transcript (sparse TSV) output. (Default: `.tsv`)
* `--out-suffix-feature-counts STR` : Suffix for the feature counts output. (Default: `.feature.counts.tsv`)
* `--out-col-id STR` : Output column name for cell ID. (Default: `cell_id`)
* `--out-col-random-key STR` : Output column name for the random key. (Default: `random_key`)
* `--out-col-x STR` / `--out-col-y STR` : Output column names for X/Y coordinates. (Defaults: `X`, `Y`)
* `--keyname-dictionary STR`, `--keyname-header-info STR`, `--keyname-n-features STR`, `--keyname-n-modalities STR`, `--keyname-n-units STR`, `--keyname-offset-data STR` : Key names used in the metadata JSON.
* `--skip-cell-tsv` : Skip generating the cell TSV output.
* `--add-xy` : Add X and Y coordinates to the output sparse TSV.
* `--seed INT` : Random seed for random key generation. (Default: 0)

## Expected Output

* `[out][out-suffix-tsv]` : The sparse TSV file (cell × feature counts).
* `[out][out-suffix-json]` : Metadata JSON (dictionary, header info, counts).
* `[out][out-suffix-cell-meta]` : Per-cell metadata table.
* `[out][out-suffix-feature-counts]` : Per-feature total counts.

## Full Usage

The full usage of `spatula pixel2sptsv` can be viewed with the `--help` option:

```
$ ./spatula pixel2sptsv --help
[./spatula pixel2sptsv] -- Convert pixel-level TSV to Sparse TSV based on assigned cell IDs

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input Options ==
   --pixel                     [STR: ]             : Input pixel-level transcript data
   --in-col-id                 [STR: cell_id]      : Column name in the input file for cell ID
   --in-col-ftr                [STR: gene]         : Column name in the input file for feature name
   --in-col-cnt                [STR: count]        : Column name in the input file for count
   --in-col-x                  [STR: X]            : Column name in the input file for x coordinate
   --in-col-y                  [STR: Y]            : Column name in the input file for y coordinate
   --idx-col-x                 [INT: 1]            : 1-based index for x coordinate column, effective only with --no-header
   --idx-col-y                 [INT: 2]            : 1-based index for y coordinate columnm, effective only with --no-header
   --idx-col-ftr               [INT: 3]            : 1-based index for feature name columnm, effective only with --no-header
   --idx-col-cnt               [INT: 4]            : 1-based index for count column, effective only with --no-header
   --idx-col-id                [INT: 5]            : 1-based index for cell ID column, effective only with --no-header
   --no-header                 [FLG: OFF]          : If set, the input file has no header line

== Input Filtering Options ==
   --ignore-ids                [STR: UNASSIGNED,-1,0,NA] : IDs to be considered as null and ignored
   --include-id-list           [STR: ]             : A file containing a list of input IDs to be included
   --exclude-id-list           [STR: ]             : A file containing a list of input IDs to be excluded
   --include-feature-list      [STR: ]             : A file containing a list of input genes to be included (feature name of IDs)
   --exclude-feature-list      [STR: ]             : A file containing a list of input genes to be excluded (feature name of IDs)
   --include-feature-regex     [STR: ]             : A regex pattern of feature/gene names to be included
   --exclude-feature-regex     [STR: ]             : A regex pattern of feature/gene names to be excluded
   --include-feature-substr    [STR: ]             : A substring of feature/gene names to be included
   --exclude-feature-substr    [STR: ]             : A substring of feature/gene names to be excluded
   --min-cell-count            [INT: 1]            : Minimum cell count to include in the output.
   --min-feature-count         [INT: 1]            : Minimum feature count to include in the output.

== Key Output Options ==
   --out                       [STR: ]             : Output prefix
   --out-suffix-json           [STR: .json]        : Suffix for the metadata JSON file
   --out-suffix-cell-meta      [STR: .cell.metadata.tsv] : Suffix for the cell metadata output file
   --out-suffix-tsv            [STR: .tsv]         : Suffix for the transcript output file
   --out-suffix-feature-counts [STR: .feature.counts.tsv] : Suffix for the feature counts output file
   --out-col-id                [STR: cell_id]      : Column name for cell ID
   --out-col-random-key        [STR: random_key]   : Column name for random key
   --out-col-x                 [STR: X]            : Column name for x coordinate
   --out-col-y                 [STR: Y]            : Column name for y coordinate
   --keyname-dictionary        [STR: dictionary]   : Key name for the dictionary in the metadata file
   --keyname-header-info       [STR: header_info]  : Key name for the header in the metadata file
   --keyname-n-features        [STR: n_features]   : Key name for the number of features in the metadata file
   --keyname-n-modalities      [STR: n_modalities] : Key name for the number of modalities in the metadata file
   --keyname-n-units           [STR: n_units]      : Key name for the number of units in the metadata file
   --keyname-offset-data       [STR: offset_data]  : Key name for the offset data in the metadata file
   --skip-cell-tsv             [FLG: OFF]          : Skip generating cell TSV output
   --add-xy                    [FLG: OFF]          : Add x and y coordinates to the output sparse TSV file

== Auxilary Output Options ==
   --seed                      [INT: 0]            : Random seed for the random key generation


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
