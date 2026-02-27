# spatula merge-sptsv

## Summary

`spatula merge-sptsv` merges multiple sparse TSV files (containing gene counts per cell/unit) into a single dataset. It handles merging of metadata (JSON), updating feature dictionaries, and can combine multiple header fields to create unique cell IDs across samples.

## Required Options

* `--list STR` : Input list file containing paths to input sparse TSV datasets. The list file can have 2 columns (`sample_id` `prefix`) or 4 columns (`sample_id` `feature_count_file` `sptsv_file` `json_file`).
* `--out STR` : Output prefix for the merged dataset.

## Additional Options

* `--combine-headers-to-id` : If enabled, multiple header fields are combined to create a unique cell ID. Useful when headers are inconsistent or need to be prefixed with sample IDs.
* `--delim STR` : Delimiter used when combining header fields. (Default: `:`)
* `--colname-combined-id STR` : Column name for the newly created combined cell ID. (Default: `cell_id`)
* `--colname-sample-id STR` : Column name for the sample ID column. (Default: `sample_id`)
* `--skip-sample-id-column` : If enabled, the sample ID column is not added to the output.
* `--min-feature-count-per-sample INT` : Minimum total count a feature must have in *each* sample to be included in the merged output. (Default: 0)
* `--colname-random-key STR` : Name of the random key column in inputs. (Default: `random_key`)
* `--exclude-random-key` : If enabled, the random key column is excluded from the output.
* `--suffix-sptsv STR` : Suffix for sparse TSV files (used when parsing 2-column list). (Default: `.tsv`)
* `--suffix-json STR` : Suffix for JSON metadata files. (Default: `.json`)
* `--suffix-feature-counts STR` : Suffix for feature count files. (Default: `.feature.counts.tsv`)

## Expected Output

The tool produces a merged dataset consisting of:
* `[out_prefix].tsv` (or specified suffix): The merged sparse TSV file containing combined data from all samples.
* `[out_prefix].json` (or specified suffix): A merged metadata file containing the updated feature dictionary and header information.
* `[out_prefix].feature.counts.tsv` (or specified suffix): A file listing the total counts of each feature in the merged dataset.

## Full Usage 

The full usage of `spatula merge-sptsv` can be viewed with the `--help` option:

```
$ ./spatula merge-sptsv --help     
[./spatula merge-sptsv] -- Merge multiple sparse TSV files into a single sparse TSV file

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --list                         [STR: ]             : Input list file containing input sparse TSV files
   --out                          [STR: ]             : Output model file to store (gene x factor) matrix

== Header Merging Options ==
   --combine-headers-to-id        [FLG: OFF]          : Combine multiple header fields to create a combined sample ID, when the headers are inconsistent across input files
   --delim                        [STR: :]            : Delimiter to use as barcode names when multiple headers fields are combined
   --colname-combined-id          [STR: cell_id]      : Column name for the combined cell ID in the output
   --colname-sample-id            [STR: sample_id]    : Column name for the sample ID in the output
   --skip-sample-id-column        [FLG: OFF]          : Skip the sample ID column in the output
   --min-feature-count-per-sample [INT: 0]            : Minimum feature count per sample to include in the output.

== Auxiliary Input/Output Options ==
   --colname-random-key           [STR: random_key]   : Column name for the random key in the output
   --keyname-dictionary           [STR: dictionary]   : Key name for the dictionary in the metadata file
   --keyname-header-info          [STR: header_info]  : Key name for the header information in the metadata file
   --keyname-n-features           [STR: n_features]   : Key name for the number of features in the metadata file
   --keyname-n-units              [STR: n_units]      : Key name for the number of units in the metadata file
   --keyname-offset-data          [STR: offset_data]  : Key name for the offset data in the metadata file
   --exclude-random-key           [FLG: OFF]          : Exclude the random key in the output
   --suffix-feature-counts        [STR: .feature.counts.tsv] : Suffix for the per-sample feature count files
   --suffix-sptsv                 [STR: .tsv]         : Suffix for the per-sample sparse TSV files
   --suffix-json                  [STR: .json]        : Suffix for the per-sample JSON metadata files


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```