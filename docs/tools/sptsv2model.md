# spatula sptsv2model

## Summary 

`spatula sptsv2model` aggregates sparse TSV data into a "model matrix" (feature x cluster) based on cluster assignments. It sums feature counts for all cells belonging to each cluster (or uses probabilistic assignments) and outputs a matrix of counts. This is useful for creating pseudobulk profiles for clusters or topics.

## Required Options

* `--out STR` : Output model file (feature x factor matrix).
* `--tsv STR` : Input TSV file containing sparse matrix data.
* `--json STR` : Input JSON metadata file.
* One of the following cluster assignment inputs is required:
    * `--clust STR` : File containing hard cluster assignments for each barcode (TSV with header).
    * `--fit STR` : File containing probabilistic cluster assignments (fit results).

## Additional Options

* `--features STR` : Input features file (name and count). Required if `--min-count` is > 0.
* `--pseudocount FLT` : Pseudocount to add to zero values. Note: The actual added value is `pseudocount / n_factors`. (Default: 1.0)
* `--min-count INT` : Minimum total count for a feature to be included in the output. (Default: 0)
* `--colname-random-key STR` : Column name for the random key in metadata. (Default: `random_key`)
* `--include-random-key` : Use the random key in barcode matching.
* `--delim STR` : Delimiter for combining header fields into barcodes. (Default: `:`)

### Metadata Keys (Advanced)
* `--keyname-dictionary STR` : Key for feature dictionary. (Default: `dictionary`)
* `--keyname-header-info STR` : Key for header info. (Default: `header_info`)
* `--keyname-n-features STR` : Key for number of features. (Default: `n_features`)
* `--keyname-n-units STR` : Key for number of units. (Default: `n_units`)
* `--keyname-offset-data STR` : Key for data offset. (Default: `offset_data`)

## Expected Output

* `[out_model]`: A TSV file where:
    * Rows correspond to features (genes).
    * Columns correspond to clusters (or factors).
    * Values are the aggregated (summed) counts for that feature in that cluster. Zero values are replaced by the pseudocount fraction.

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