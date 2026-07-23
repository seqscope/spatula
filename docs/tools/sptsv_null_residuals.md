# spatula sptsv-null-residuals

## Summary

`spatula sptsv-null-residuals` quantifies how far a sparse TSV (FICTURE2 format) deviates from a null model in which gene identity and unit identity are independent. Under that null model, the expected count of gene $g$ in unit $u$ is

```
expected[g,u] = marginal[g] * marginal[u] / total
```

where `marginal[g]` is the total count of gene $g$ across all units (taken from the `.feature.stats.tsv` file), `marginal[u]` is the total count of unit $u$ across all genes (the total transcript count of the unit, i.e. the 5th column of the `.txt` file when `offset_data` is 3), and `total` is the grand total of all counts.

The tool reports the L1 deviation from that expectation, summed per unit and per feature:

```
residual[u] = sum over all genes of |x[g,u] - expected[g,u]|
residual[g] = sum over all units of |x[g,u] - expected[g,u]|
```

Both sums run over **all** gene/unit pairs, including the zero entries that the sparse format does not store. Because `sum_g expected[g,u] = marginal[u]` and `sum_u expected[g,u] = marginal[g]`, the contribution of the zero entries has a closed form, so the residuals are computed exactly in a single streaming pass over the sparse TSV.

The residual of a unit is bounded above by `2 * marginal[u]` (and likewise for a feature), so the reported `frac_residual` is always between 0 (the unit/feature is exactly what the null model predicts) and 1 (maximal deviation).

## Required Options

* `--sptsv STR` : Input sparse TSV prefix. The tool reads `[prefix].txt`, `[prefix].json`, and `[prefix].feature.stats.tsv`.
* `--out STR` : Output prefix for the per-unit and per-feature residual files.

## Additional Options

* `--min-feature-count INT` : Minimum total count for a feature to be included in the null model. Excluded features are dropped from both outputs and from the marginals/total, so unit marginals are recomputed from the retained features only. (Default: 0)
* `--exclude-random-key` : Do not carry the random key column over to the per-unit output.
* `--colname-random-key STR` : Column name of the random key in the metadata. (Default: `random_key`)
* `--suffix-sptsv-tsv STR` : Suffix of the input sparse TSV file. (Default: `.txt`)
* `--suffix-sptsv-json STR` : Suffix of the input JSON metadata file. (Default: `.json`)
* `--suffix-feature-stats STR` : Suffix of the input feature stats file. (Default: `.feature.stats.tsv`)
* `--suffix-unit-stats STR` : Suffix of the output per-unit file. Use a `.gz` suffix to write a compressed file. (Default: `.unit_stats_null.tsv`)
* `--suffix-feature-residuals STR` : Suffix of the output per-feature file. Use a `.gz` suffix to write a compressed file. (Default: `.feature_residuals_null.tsv`)

### Metadata Keys (Advanced)
* `--keyname-dictionary STR` : Key for feature dictionary. (Default: `dictionary`)
* `--keyname-header-info STR` : Key for header info. (Default: `header_info`)
* `--keyname-n-features STR` : Key for number of features. (Default: `n_features`)
* `--keyname-n-units STR` : Key for number of units. (Default: `n_units`)
* `--keyname-offset-data STR` : Key for data offset. (Default: `offset_data`)

## Expected Output

* `[out].unit_stats_null.tsv` : One row per unit, with a header line. The leading columns are the header columns of the input sparse TSV (as listed in `header_info`, e.g. `random_key`, `X`, `Y`), followed by:
    * `n_features` : Number of distinct retained features observed in the unit.
    * `n_umis` : `marginal[u]`, the total count of the unit over the retained features.
    * `l1_residual` : `sum over all genes of |x[g,u] - expected[g,u]|`.
    * `frac_residual` : `l1_residual / (2 * n_umis)`, in [0, 1].
* `[out].feature_residuals_null.tsv` : One row per retained feature, with a header line:
    * `feature` : Feature (gene) name.
    * `feature_idx` : Feature index used in the sparse TSV.
    * `n_units` : Number of units with a non-zero count for the feature.
    * `n_umis` : `marginal[g]`, the total count of the feature.
    * `l1_residual` : `sum over all units of |x[g,u] - expected[g,u]|`.
    * `frac_residual` : `l1_residual / (2 * n_umis)`, in [0, 1].

The sum of `l1_residual` over all units equals the sum over all features, because both are the total L1 deviation of the matrix from its null expectation.

The tool also cross-checks the inputs and emits warnings (without stopping) if the per-unit total does not match the 5th column of the `.txt` file, if the observed per-feature totals do not match the `.feature.stats.tsv` file, or if the unit and feature marginals do not sum to the same grand total.

## Full Usage

The full usage of `spatula sptsv-null-residuals` can be viewed with the `--help` option:

```
$ ./spatula sptsv-null-residuals --help
[./spatula sptsv-null-residuals] -- Compute per-unit and per-feature L1 residuals from the null (independence) model of a Sparse TSV

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --sptsv                    [STR: ]             : Input sparse TSV prefix
   --out                      [STR: ]             : Output prefix
   --min-feature-count        [INT: 0]            : Minimum total count for a feature to be included in the null model

== Auxiliary Input/Output Options ==
   --colname-random-key       [STR: random_key]   : Column name for the random key in the output
   --keyname-dictionary       [STR: dictionary]   : Key name for the dictionary in the metadata file
   --keyname-header-info      [STR: header_info]  : Key name for the header information in the metadata file
   --keyname-n-features       [STR: n_features]   : Key name for the number of features in the metadata file
   --keyname-n-units          [STR: n_units]      : Key name for the number of units in the metadata file
   --keyname-offset-data      [STR: offset_data]  : Key name for the offset data in the metadata file
   --exclude-random-key       [FLG: OFF]          : Exclude the random key column in the per-unit output
   --suffix-feature-residuals [STR: .feature_residuals_null.tsv] : Suffix for the output per-feature residual file
   --suffix-unit-stats        [STR: .unit_stats_null.tsv] : Suffix for the output per-unit residual file
   --suffix-sptsv-json        [STR: .json]        : Suffix for the input JSON metadata file
   --suffix-feature-stats     [STR: .feature.stats.tsv] : Suffix for the input feature stats file
   --suffix-sptsv-tsv         [STR: .txt]         : Suffix for the input sparse TSV file


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
