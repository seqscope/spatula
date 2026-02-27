# spatula merge-pseudobulk

## Summary 

`spatula merge-pseudobulk` merges multiple pseudobulk count matrices (e.g. outputs from LDA models or other aggregation tools) into a single matrix. It verifies consistency of row names (genes/features) and column names (topics/samples) before summing the counts.

## Required Options

* `--out STR` : Output filename for the merged pseudobulk matrix.
* One of the following input options must be specified:
    * `--tsv STR` : Input pseudobulk TSV file. Can be specified multiple times to merge multiple files.
    * `--list STR` : File containing a list of input pseudobulk TSV files to merge (one per line).

## Additional Options

N/A

## Expected Output

* A single TSV file containing the merged pseudobulk matrix, where counts from input matrices are summed for matching features and samples. The output format matches the input pseudobulk TSV format (typically genes as rows, topics/clusters as columns).

## Full Usage 

The full usage of `spatula merge-pseudobulk` can be viewed with the `--help` option:

```
$ ./spatula merge-pseudobulk --help 
[./spatula merge-pseudobulk] -- Merge multiple pseudobulk matrices into one

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input/Output options ==
   --list [STR: ]             : File containing list of pseudobulk tsv files to merge
   --tsv  [V_STR: ]           : Input pseudobulk tsv files to merge (can be specified multiple times)
   --out  [STR: ]             : Output filename


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```