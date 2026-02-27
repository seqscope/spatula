# spatula sptsv2mex

## Summary 

`spatula sptsv2mex` converts the Sparse TSV format (used by FICTURE) into the 10x Gene Expression Marker (MEX) format (Matrix Market format). This conversion allows downstream analysis with tools that expect standard 10x output formats.

## Required Options

* `--tsv STR` : Input TSV file containing the sparse matrix data (e.g., as output by `mex2sptsv` or `merge-sptsv`).
* `--json STR` : Input JSON metadata file containing header information, feature dictionary, etc.
* `--out-dir STR` : Output directory where the MEX files (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`) will be saved.

## Additional Options

* `--bcd STR` : Filename for the output barcode file. (Default: `barcodes.tsv.gz`)
* `--ftr STR` : Filename for the output feature file. (Default: `features.tsv.gz`)
* `--mtx STR` : Filename for the output matrix file. (Default: `matrix.mtx.gz`)
* `--colname-random-key STR` : Column name for the random key in the JSON metadata. (Default: `random_key`)
* `--include-random-key` : If enabled, the random key column is included in the generated barcode names.
* `--delim STR` : Delimiter used when combining multiple header fields into a single barcode string. (Default: `:`)
* `--feature-type STR` : Feature type label to use in the output `features.tsv.gz` file (e.g., "Gene Expression"). (Default: `Gene Expression`)

### Metadata Keys (Advanced)
These options allow customization if the input JSON uses non-standard keys.
* `--keyname-dictionary STR` : Key for the feature dictionary. (Default: `dictionary`)
* `--keyname-header-info STR` : Key for the header information. (Default: `header_info`)
* `--keyname-n-features STR` : Key for the number of features. (Default: `n_features`)
* `--keyname-n-units STR` : Key for the number of units/cells. (Default: `n_units`)
* `--keyname-offset-data STR` : Key for the data offset. (Default: `offset_data`)

## Expected Output

The tool creates the specified output directory and populates it with three files standard to the 10x MEX format:
1.  **`barcodes.tsv.gz`**: List of barcodes (cell IDs).
2.  **`features.tsv.gz`**: List of features (genes), including ID, name, and type.
3.  **`matrix.mtx.gz`**: The sparse expression matrix in Matrix Market format.

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