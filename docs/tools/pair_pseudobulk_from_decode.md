# spatula pair-pseudobulk-from-decode

## Summary

`spatula pair-pseudobulk-from-decode` writes a *pairwise* pseudobulk matrix from pixel-level decode output. For each transcript it reads two factor assignments (`--colname-factor1` and `--colname-factor2`) and aggregates feature counts across each pair of factors, producing both a gene-by-factor-pair table and a factor-pair summary table. This is a pairwise generalization of [pseudobulk-from-decode](pseudobulk_from_decode.md), useful for examining co-occurrence of factors.

## Required Options

* `--tsv STR` : Input pixel-level TSV file (use `/dev/stdin` to read from standard input).
* `--out STR` : Output prefix.
* `--colname-factor1 STR` : Column name for the first factor.
* `--colname-factor2 STR` : Column name for the second factor.

## Additional Options

* `--colname-feature STR` : Column name for the feature. (Default: `gene`)
* `--colname-count STR` : Column name for the count. (Default: `count`)
* `--missing-value-str STR` : String representing a missing value. (Default: `NA`)
* `--out-g-f1-f2-suffix STR` : Suffix for the gene × (factor1, factor2) output file. (Default: `.g_f1_f2.tsv.gz`)
* `--out-f1-f2-suffix STR` : Suffix for the (factor1, factor2) summary output file. (Default: `.f1_f2.tsv.gz`)

## Expected Output

* `[out][out-g-f1-f2-suffix]` : Pseudobulk counts per gene for each pair of factors.
* `[out][out-f1-f2-suffix]` : Summary counts for each pair of factors.

## Full Usage

The full usage of `spatula pair-pseudobulk-from-decode` can be viewed with the `--help` option:

```
$ ./spatula pair-pseudobulk-from-decode --help
[./spatula pair-pseudobulk-from-decode] -- Write pairwise pseudobulk matrix from pixel-level decode output

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --tsv                [STR: ]             : tsv file to draw the x-y coordinates. /dev/stdin for stdin
   --colname-feature    [STR: gene]         : Column name for the feature
   --colname-count      [STR: count]        : Column name for the count
   --colname-factor1    [STR: ]             : Column name for the first factor
   --colname-factor2    [STR: ]             : Column name for the second factor
   --missing-value-str  [STR: NA]           : Missing value string

== Output Options ==
   --out                [STR: ]             : Output prefix
   --out-g-f1-f2-suffix [STR: .g_f1_f2.tsv.gz] : Suffix for the output TSV file (default: .g_f1_f2.tsv.gz)
   --out-f1-f2-suffix   [STR: .f1_f2.tsv.gz] : Suffix for the output TSV file (default: .f1_f2.tsv.gz)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
