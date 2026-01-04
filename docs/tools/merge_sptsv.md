# spatula merge-sptsv

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

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