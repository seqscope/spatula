# spatula sptsv2model

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

## Full Usage 

The full usage of `spatula sptsv2model` can be viewed with the `--help` option:

```
$ ./spatula sptsv2model --help           
[./spatula sptsv2model] -- Create model matrix from Sparse TSV format in FICTURE2 with cluster assignment

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input Options ==
   --tsv                 [STR: ]             : Input TSV file containing the sparse matrix data
   --json                [STR: ]             : Input JSON metadata file containing the header information
   --features            [STR: ]             : Input features file containing the feature names and counts
   --clust               [STR: ]             : Input cluster file containing the cluster assignments for each barcode
   --fit                 [STR: ]             : Input fit results file containing the probabilistic cluster assignment for each barcode

== Key Output Options ==
   --out                 [STR: ]             : Output model file to store (gene x factor) matrix

== Auxiliary Input/Output Options ==
   --pseudocount         [FLT: 1.00]         : Pseudocount to use when normalizing the feature counts (default: 1.0)
   --min-count           [INT: 0]            : Minimum feature count to include in the output.
   --colname-random-key  [STR: random_key]   : Column name for the random key in the output
   --include-random-key  [FLG: OFF]          : Include the random key in the output
   --keyname-dictionary  [STR: dictionary]   : Key name for the dictionary in the metadata file
   --keyname-header-info [STR: header_info]  : Key name for the header information in the metadata file
   --keyname-n-features  [STR: n_features]   : Key name for the number of features in the metadata file
   --keyname-n-units     [STR: n_units]      : Key name for the number of units in the metadata file
   --keyname-offset-data [STR: offset_data]  : Key name for the offset data in the metadata file
   --delim               [STR: :]            : Delimiter to use as barcode names when multiple headers fields are combined


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```