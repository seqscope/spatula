# spatula pseudobulk-from-decode

## Summary 

`spatula pseudobulk-from-decode` generates a pseudobulk expression matrix from pixel-level decode output. It aggregates counts for each feature across spatial locations defined by K1, K2, K3 factors.

An example command is given below:

```bash
spatula pseudobulk-from-decode --tsv /path/to/pixel.decode.tsv.gz --colname-count count --colname-feature gene --colname-K1 t24-f48-p24-a6_K1 --colname-P1 "" --colname-K2 "" --colname-P2 "" --colname-K3 "" --colname-P3 "" --write-cell-tsv --out /path/to/output.prefix.pseudobulk.tsv.gz --out-cell-tsv /path/to/output.prefix.cells.tsv.gz 
```

## Required options

* `--tsv`: Input TSV file containing pixel-level decode output.
* `--out` : Output file name for the pseudobulk expression matrix (TSV).

## Options for writing per-cell results matrix

* `--write-cell-tsv` : If set, writes a per-cell factor assignment distribution TSV file.
* `--out-cell-tsv` : Output file name for the per-cell results TSV file
* `--colname-cell` : Column name for the cell ID in the per-cell results TSV file

## Other Additional Options

* `--colname-feature`: Column name for the feature (default: `feature`).
* `--colname-count`: Column name for the count (default: `ct`).
* `--colname-K1`: Column name for the K1 value (default: `K1`).
* `--colname-K2`: Column name for the K2 value (default: `K2`).
* `--colname-K3`: Column name for the K3 value (default: `K3`).
* `--colname-P1`: Column name for the P1 value (default: `P1`).
* `--colname-P2`: Column name for the P2 value (default: `P2`).
* `--colname-P3`: Column name for the P3 value (default: `P3`).
* `--n-factors`: Force the total number of factors. Only works with integer factors.

## Expected Output

The output pseudobulk expression matrix TSV file contains rows for each feature and columns for each factor (combination of K1, K2, K3). Each entry represents the aggregated count for that feature in the corresponding factor. The first column contains the feature (gene) names, and the subsequent columns contain the counts for each factor.

The output per-cell results TSV file (if `--write-cell-tsv` is set) contains the distribution of each factor (sum to 1) for each cell. The first column contains the cell IDs, and the subsequent columns contain the probabilities for each factor.

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