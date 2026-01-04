# spatula pseudobulk-from-decode

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

## Full Usage 

The full usage of `spatula pseudobulk-from-decode` can be viewed with the `--help` option:

```
$ ./spatula pseudobulk-from-decode --help
[./spatula pseudobulk-from-decode] -- Write pseudobulk matrix from pixel-level decode output

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --tsv             [STR: ]             : tsv file to draw the x-y coordinates. /dev/stdin for stdin
   --colname-feature [STR: feature]      : Column name for the feature
   --colname-count   [STR: ct]           : Column name for the count
   --colname-K1      [STR: K1]           : Column name for the K1 value
   --colname-K2      [STR: K2]           : Column name for the K2 value
   --colname-K3      [STR: K3]           : Column name for the K3 value
   --colname-P1      [STR: P1]           : Column name for the P1 value
   --colname-P2      [STR: P2]           : Column name for the P2 value
   --colname-P3      [STR: P3]           : Column name for the P3 value

== Output Options ==
   --out             [STR: ]             : Output file name
   --n-factors       [INT: 0]            : Force the total number of factors. Only works with integer factors. Encoded as 0, 1, 2, ..., (n_factors-1)

== Auxiliary Options to generate per-cell results matrix ==
   --write-cell-tsv  [FLG: OFF]          : Write per-cell results TSV file
   --out-cell-tsv    [STR: ]             : Output file name for per-cell results
   --colname-cell    [STR: cell_id]      : Column name for the cell ID


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run

```