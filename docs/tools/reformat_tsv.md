# spatula reformat-tsv

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

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