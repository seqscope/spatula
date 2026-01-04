# spatula sptsv2mex

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

## Full Usage 

The full usage of `spatula sptsv2mex` can be viewed with the `--help` option:

```
$ ./spatula sptsv2mex --help  
[./spatula sptsv2mex] -- Convert Sparse TSV format in FICTURE2 to 10x Market Exchange (MEX) formatted DGE

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input Options ==
   --tsv                 [STR: ]             : Input TSV file containing the sparse matrix data
   --json                [STR: ]             : Input JSON metadata file containing the header information

== Key Output Options ==
   --out-dir             [STR: ]             : Output directory for the MEX files
   --bcd                 [STR: barcodes.tsv.gz] : Barcode file name
   --ftr                 [STR: features.tsv.gz] : Feature file name
   --mtx                 [STR: matrix.mtx.gz] : Matrix file name

== Auxiliary Input/Output Options ==
   --colname-random-key  [STR: random_key]   : Column name for the random key in the output
   --include-random-key  [FLG: OFF]          : Include the random key in the output
   --keyname-dictionary  [STR: dictionary]   : Key name for the dictionary in the metadata file
   --keyname-header-info [STR: header_info]  : Key name for the header information in the metadata file
   --keyname-n-features  [STR: n_features]   : Key name for the number of features in the metadata file
   --keyname-n-units     [STR: n_units]      : Key name for the number of units in the metadata file
   --keyname-offset-data [STR: offset_data]  : Key name for the offset data in the metadata file
   --delim               [STR: :]            : Delimiter to use as barcode names when multiple headers fields are combined
   --feature-type        [STR: Gene Expression] : Feature type to use in the output files (default: Gene Expression)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```