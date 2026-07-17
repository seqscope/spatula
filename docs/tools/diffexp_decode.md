# spatula diffexp-decode

## Summary

`spatula diffexp-decode` performs a differential expression test and creates pseudobulk data directly from pixel-level decoding output (e.g. from FICTURE). For each pixel, the top factor assignments (`K1`/`K2`/`K3`) and their posterior probabilities (`P1`/`P2`/`P3`) are used to aggregate feature counts into per-factor pseudobulk profiles, and a chi-square-based differential expression test is applied to identify features enriched in each factor.

## Required Options

* `--tsv STR` : Input pixel-level TSV file (use `/dev/stdin` to read from standard input).
* `--out STR` : Output file prefix.

## Additional Options

### Input Column Names
* `--colname-feature STR` : Column name for the feature. (Default: `feature`)
* `--colname-count STR` : Column name for the count. (Default: `ct`)
* `--colname-K1 STR` / `--colname-K2 STR` / `--colname-K3 STR` : Column names for the top-1/2/3 factor assignments. (Defaults: `K1`, `K2`, `K3`)
* `--colname-P1 STR` / `--colname-P2 STR` / `--colname-P3 STR` : Column names for the top-1/2/3 posterior probabilities. (Defaults: `P1`, `P2`, `P3`)

### Test Settings
* `--max-pval FLT` : Maximum p-value for a feature to be reported as differentially expressed. (Default: `1.0e-03`)
* `--min-fc FLT` : Minimum fold change for the differential expression test. (Default: `1.50`)
* `--pseudocount FLT` : Pseudocount added in the differential expression test. (Default: `0.50`)

### Output Suffixes
* `--suffix_de STR` : Suffix for the differential expression output. (Default: `.bulk_chisq.tsv`)
* `--suffix_post STR` : Suffix for the posterior count (pseudobulk) output. (Default: `.posterior.count.tsv.gz`)

## Expected Output

* `[out][suffix_de]` : Differential expression results (per-factor chi-square test statistics, fold changes, and p-values).
* `[out][suffix_post]` : Posterior (pseudobulk) count matrix aggregated per factor.

## Full Usage

The full usage of `spatula diffexp-decode` can be viewed with the `--help` option:

```
$ ./spatula diffexp-decode --help
[./spatula diffexp-decode] -- Perform differential expression test and create pseudobulk data from pixel-level decoding output

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

== Settings ==
   --max-pval        [FLT: 1.0e-03]      : Max p-value for the differential expression test
   --min-fc          [FLT: 1.50]         : Min fold change for the differential expression test
   --pseudocount     [FLT: 0.50]         : Pseudocount for the differential expression test

== Output Options ==
   --out             [STR: ]             : Output file name
   --suffix_de       [STR: .bulk_chisq.tsv] : Suffix for differential expression output
   --suffix_post     [STR: .posterior.count.tsv.gz] : Suffix for posterior count output


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
