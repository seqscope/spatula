# spatula mex2sptsv

## Summary 

`spatula mex2sptsv` converts 10x Gene Expression Marker (MEX) format (Matrix Market format) into the Sparse TSV format used by FICTURE. It takes a directory containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` (or specified filenames) and outputs a sparse TSV file along with metadata and feature count files. It supports filtering barcodes and features during conversion.

## Required Options

* `--out STR` : Output prefix for the generated files.

## Additional Options

### Input Options
* `--in-dir STR` : Input directory containing the MEX files.
* `--bcd STR` : Filename for barcodes. (Default: `barcodes.tsv.gz`)
* `--ftr STR` : Filename for features. (Default: `features.tsv.gz`)
* `--mtx STR` : Filename for the matrix. (Default: `matrix.mtx.gz`)
* `--icol-mtx INT` : 1-based column index in the matrix file to use as the count. (Default: 3)

### Filtering Options
* `--include-feature-list STR` : File with list of feature IDs/names to include.
* `--exclude-feature-list STR` : File with list of feature IDs/names to exclude.
* `--include-feature-regex STR` : Regex pattern for features to include.
* `--exclude-feature-regex STR` : Regex pattern for features to exclude.
* `--include-feature-substr STR` : Substring for features to include.
* `--exclude-feature-substr STR` : Substring for features to exclude.
* `--include-barcode-list STR` : File with list of barcode IDs to include.
* `--exclude-barcode-list STR` : File with list of barcode IDs to exclude.
* `--min-feature-count INT` : Minimum total count for a feature to be included. (Default: 0)

### Output Configuration
* `--out-suffix-meta STR` : Suffix for metadata JSON. (Default: `.json`)
* `--out-suffix-tsv STR` : Suffix for sparse TSV. (Default: `.tsv`)
* `--out-suffix-feature-counts STR` : Suffix for feature counts. (Default: `.feature.counts.tsv`)
* `--colname-barcode-name STR` : Column name for barcode ID. (Default: `cell_id`)
* `--colname-random-key STR` : Column name for random key. (Default: `random_key`)
* `--seed INT` : Seed for random key generation. (Default: 0, which uses time)
* `--use-gene-id` : Use gene IDs instead of gene names in output.
* `--error-on-duplicate-gene-name` : Fail if duplicate gene names are found. Otherwise, duplicates are renamed.

## Expected Output

* `[out_prefix].tsv` (or specified suffix): The converted sparse TSV file.
* `[out_prefix].json` (or specified suffix): Metadata file including feature dictionary and header info.
* `[out_prefix].feature.counts.tsv` (or specified suffix): File containing total counts for each feature.

## Full Usage 

The full usage of `spatula mex2sptsv` can be viewed with the `--help` option:

```
$ ./spatula mex2sptsv --help  
[./spatula mex2sptsv] -- Convert 10x Market Exchange (MEX) formatted DGE to Sparse TSV format in FICTURE2

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input Options ==
   --in-dir                       [STR: ]             : Input directory containing the MEX files
   --bcd                          [STR: barcodes.tsv.gz] : Barcode file name
   --ftr                          [STR: features.tsv.gz] : Feature file name
   --mtx                          [STR: matrix.mtx.gz] : Matrix file name
   --icol-mtx                     [INT: 3]            : 1-based column index in the matrix file to use as the count

== Input Filtering Options ==
   --include-feature-list         [STR: ]             : A file containing a list of input genes to be included (feature name of IDs)
   --exclude-feature-list         [STR: ]             : A file containing a list of input genes to be excluded (feature name of IDs)
   --include-feature-regex        [STR: ]             : A regex pattern of feature/gene names to be included
   --exclude-feature-regex        [STR: ]             : A regex pattern of feature/gene names to be excluded
   --include-feature-substr       [STR: ]             : A substring of feature/gene names to be included
   --exclude-feature-substr       [STR: ]             : A substring of feature/gene names to be excluded
   --include-barcode-list         [STR: ]             : A file containing a list of input barcode IDs to be included (feature name of IDs)
   --exclude-barcode-list         [STR: ]             : A file containing a list of input barcode IDs to be excluded (feature name of IDs)
   --min-feature-count            [INT: 0]            : Minimum feature count to include in the output. Requires reading the input matrix twice

== Key Output Options ==
   --out                          [STR: ]             : Output prefix
   --out-suffix-meta              [STR: .json]        : Suffix for the metadata output file
   --out-suffix-tsv               [STR: .tsv]         : Suffix for the transcript output file
   --out-suffix-feature-counts    [STR: .feature.counts.tsv] : Suffix for the feature counts output file
   --colname-barcode-name         [STR: cell_id]      : Column name for barcode name
   --colname-random-key           [STR: random_key]   : Column name for random key
   --keyname-dictionary           [STR: dictionary]   : Key name for the dictionary in the metadata file
   --keyname-header-info          [STR: header_info]  : Key name for the header in the metadata file
   --keyname-n-features           [STR: n_features]   : Key name for the number of features in the metadata file
   --keyname-n-modalities         [STR: n_modalities] : Key name for the number of modalities in the metadata file
   --keyname-n-units              [STR: n_units]      : Key name for the number of units in the metadata file
   --keyname-offset-data          [STR: offset_data]  : Key name for the offset data in the metadata file

== Auxilary Output Options ==
   --seed                         [INT: 0]            : Random seed for the random key generation
   --use-gene-id                  [FLG: OFF]          : Use gene ID instead of gene name in the output
   --error-on-duplicate-gene-name [FLG: OFF]          : Report error on duplicate gene name in the feature file. If not set, the duplicate gene names will be renamed by combining with the gene IDs


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```