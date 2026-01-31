# spatula reformat-tsv

## Summary 

`spatula reformat-tsv` is a flexible tool for reformatting TSV/CSV files. It allows selecting specific columns, renaming them, adding constant columns, scaling numeric values, and changing delimiters. It can also compute min/max statistics for coordinates and aggregate feature counts.

## Required Options

* `--in STR` : Input CSV/TSV file.
* `--out STR` : Output filename prefix.

## Additional Options

### Input/Output Formatting
* `--in-delim STR` : Delimiter character for the input file. (Default: `\t`)
* `--out-delim STR` : Delimiter character for the output file. (Default: `\t`)
* `--unquote` : Remove quotes from string values in the input file.
* `--skip-lines INT` : Number of lines to skip at the beginning of the input file. (Default: 0)

### Column Selection and Modification
* `--colnames STR` : Comma-separated list of columns to include in the output.
    * Use `col` to include a column.
    * Use `old:new` to rename a column.
    * Use `:name:value` to add a new column `name` with a constant `value`.
* `--include-rest` : If enabled, include all remaining columns from the input that were not specified in `--colnames`.
* `--scale V_STR` : Scale numeric columns. Format: `colname:factor[:format]`.
    * `factor`: Multiplier.
    * `format`: Output format (e.g., `2f` for `%.2f`, `3e` for `%.3e`). Default is `f`.
* `--add-count` : (Flag present in help but logic seems missing/unused in main loop? Code has `add_count` param but not used in loop). *Note based on code: `add_count` is parsed but seemingly not implemented in the loop.*

### Auxiliary Outputs
* `--write-minmax` : Write a file containing min/max values for X and Y coordinates.
    * Requires `--colname-x` and `--colname-y` to identify coordinate columns.
* `--write-features` : Write a file containing unique feature counts.
    * Requires `--colname-feature` and optionally `--colname-count`.
* `--colname-x STR` : Name of X column (for minmax). (Default: `X`)
* `--colname-y STR` : Name of Y column (for minmax). (Default: `Y`)
* `--colname-feature STR` : Name of feature column (for feature counts). (Default: `gene`)
* `--colname-count STR` : Name of count column (for feature counts). (Default: `count`)
* `--suffix-tsv STR` : Suffix for the main output file. (Default: `.tsv.gz`)
* `--suffix-minmax STR` : Suffix for the minmax output file. (Default: `.minmax.tsv`)
* `--suffix-features STR` : Suffix for the feature counts output file. (Default: `.features.tsv.gz`)

## Expected Output

* `[out_prefix][suffix-tsv]`: The reformatted TSV/CSV file.
* `[out_prefix][suffix-minmax]` (Optional): A file with `xmin`, `xmax`, `ymin`, `ymax` values.
* `[out_prefix][suffix-features]` (Optional): A file with feature names and their total counts.

## Full Usage 

The full usage of `spatula reformat-tsv` can be viewed with the `--help` option:

```
$ ./spatula reformat-tsv --help
[./spatula reformat-tsv] -- Reformat TSV/CSV files by selecting or reordering columns

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --in              [STR: ]             : Input CSV/TSV file
   --in-delim        [STR: 	]            : Input delimiter
   --colnames        [STR: ]             : Comma-separated column names to include in the output. To rename columns, use 'original_name:new_name' format. To add constants, use ':name:value' format.
   --scale           [V_STR: ]           : Scale the columns by multiplying with the given value. Use 'colname:scale:format' format. use '0f' for integer, or use '3f', '3e', '5g', etc. Use the original column name before renaming
   --include-rest    [FLG: OFF]          : Include all columns not specified in --colnames at the end
   --unquote         [FLG: OFF]          : Unquote the string values in the input file
   --skip-lines      [INT: 0]            : Number of lines to skip at the beginning of the input file

== Output Options ==
   --out             [STR: ]             : Output CSV/TSV file prefix
   --out-delim       [STR: 	]            : Output delimiter
   --suffix-tsv      [STR: .tsv.gz]      : Output suffix for TSV files
   --suffix-minmax   [STR: .minmax.tsv]  : Output suffix for minmax files
   --suffix-features [STR: .features.tsv.gz] : Output suffix for features files
   --colname-x       [STR: X]            : Column name for x-axis
   --colname-y       [STR: Y]            : Column name for y-axis
   --colname-feature [STR: gene]         : Column name for feature
   --colname-count   [STR: count]        : Column name for count
   --write-minmax    [FLG: OFF]          : Write minmax file
   --write-features  [FLG: OFF]          : Write features file
   --add-count       [FLG: OFF]          : Add count column to the output file, assigning 1 to all rows


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```