# spatula diffexp-model-matrix

## Summary 

`spatula diffexp-model-matrix` performs differential expression testing given model matrix TSV files. When one model matrix TSV file is provided, it performs per-factor differential expression tests for each gene. When two model matrix TSV files are provided, it performs both bulk differential expression test as well as per-factor conditional differential expression tests for each gene. Only simple chi-square test is currently implemented.

An example usage of the tool is as follows:

```sh
## Usage for performing differential expression test between two model matrices
spatula diffexp-model-matrix --tsv1 /path/to/input.model1.matrix.tsv.gz --tsv2 /path/to/input.model2.matrix.tsv.gz \
                    --out /path/to/output.prefix --min-fc 1.5 --max-pval 1e-3 --min-count 10 --pseudocount 0.5
```

See below for a more detailed usage description.

## Required options

* `--tsv1 STR` : First TSV file containing the model or pseudobulk matrix where the first column is the feature (gene) name followed by per-topic counts. If only one TSV file is provided, per-factor differential expression tests will be performed.
* `--out STR` : Output file prefix for storing the differential expression test results.

## Additional Options

* `--tsv2 STR` : Second TSV file containing the model or pseudobulk matrix. If provided, the two matrices will be compared for differential expression testing.
* `--max-pval FLT [1.0e-03]` : Maximum p-value threshold for reporting differential expression test results. Default is 1.0e-03.
* `--min-fc FLT [1.50]` : Minimum fold change threshold for reporting differential expression test results. Default is 1.50.
* `--min-count FLT [10.00]` : Minimum observed count for the feature to be considered to report differential expression test results. Default is 10.00.
* `--pseudocount FLT [0.50]` : Pseudocount to add to the counts to avoid zero counts. Default is 0.50.
* `--test-pairwise FLG [OFF]` : If set, perform pairwise test between all possible pairs of factors (only applicable when one TSV file is provided). Default is OFF.
* `--ignore-mismatch FLG [OFF]` : If set, ignore mismatching factors between the two TSV files. Only overlapping factors will be used for the differential expression test. If not set, an error will be raised when mismatching factors are found. Default is OFF.

## Expected Output

### With one input TSV file

The following output files will be generated based on the provided options:
* `{outprefix}.de.marginal.tsv.gz` : Marginal differential expression test results per feature (gene)
  * When single TSV file is provided, this file contains per-factor differential expression test results, comparing each factor against all other factors combined for each gene. The 2x2 contigency table is constructed for each gene and factor as follows:

|               | Factor i | Other Factors |
|---------------|----------|---------------|
| Gene g Count  |    a     |      b        |
| Other Genes   |    c     |      d        |

* `{outprefix}.de.pairwise.tsv.gz` : When `--test-pairwise` option is ON, pairwise differential expression test results per feature (gene) between all possible pairs of factors will be saved in this file. The 2x2 contigency table is constructed for each factor pair `(i, j)` per gene as follows:

|               | Factor i |  Factor j  |
|---------------|----------|------------|
| Gene g Count  |    a     |      b     |
| Other Genes   |    c     |      d     |

### With two input TSV files

The following output files will be generated based on the provided options:

* `{outprefix}.de.bulk.feature.tsv.gz` : Bulk differential expression test results per feature (gene) between the two TSV files, summing up the counts across all factors. The 2x2 contigency table is constructed for each gene as follows:

|               | TSV1 | TSV2 |
|---------------|------|------|
| Gene g Count  |  a   |  b   |
| Other Genes   |  c   |  d   |

where `a` is the count of gene `g` in factor `i`, `b` is the count of gene `g` in all other factors, `c` is the count of all other genes in factor `i`, and `d` is the count of all other genes in all other factors.

* `{outprefix}.de.bulk.factor.tsv.gz` : Bulk differential expression test results per factor between the two TSV files, summing up the counts across all genes. The 2x2 contigency table is constructed for each gene as follows:

|                 | TSV1 | TSV2 |
|-----------------|------|------|
| Factor f Count  |  a   |  b   |
| Other Factors   |  c   |  d   |

where `a` is the count of factor `f` in TSV1, `b` is the count of factor `f` in TSV2, `c` is the count of all other factors in TSV1, and `d` is the count of all other factors in TSV2. This test identifies factors (e.g. cell types) that have differential abundances between the two TSV files.


* `{outprefix}.de.conditional.feature.tsv.gz` : Conditional differential expression test results per feature (gene) between the two TSV files for each factor, controlling for the other factors. The 2x2 contigency table is constructed for each gene and factor as follows:

|                           | TSV1 | TSV2 |
|---------------------------|------|------|
| Gene g Count for Factor f |  a   |  b   |
| Other genes for Factor f  |  c   |  d   |

The output contains the following common columns:
* `Chi2` : Chi-square statistic
* `pval` : P-value from the chi-square test (can be zero if underflow occurs)
* `log10pval` : -log10 (P-value) 
* `FoldChange` : Fold change between the foreground vs background (with one TSV file).
* `log2FC` : log2 Fold change between TSV1 and TSV2 (with two TSV files).


## Full Usage 

The full usage of `spatula diffexp-model-matrix` can be viewed with the `--help` option:

```
$ ./spatula diffexp-model-matrix --help
[./spatula diffexp-model-matrix] -- Perform differential expression test from model matrix

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input/Output options ==
   --tsv1            [STR: ]             : First tsv file containing the model or pseudobulk matrix
   --tsv2            [STR: ]             : (Optional) Second tsv file containing the model or pseudobulk matrix
   --out             [STR: ]             : Output file prefix

== Settings ==
   --max-pval        [FLT: 1.0e-03]      : Max p-value for the differential expression test
   --min-fc          [FLT: 1.50]         : Min fold change for the differential expression test
   --min-count       [FLT: 10.00]        : Minimum observed count for the feature to be considered
   --pseudocount     [FLT: 0.50]         : Pseudocount to add to the counts
   --test-pairwise   [FLG: OFF]          : Perform pairwise test (1 sample test only)
   --ignore-mismatch [FLG: OFF]          : Ignore mismatching factors between tsv1 and tsv2. If set, only overlapping factors will be used for the test


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run

```