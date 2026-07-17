# spatula assign-feature2bin

## Summary

`spatula assign-feature2bin` assigns each feature (gene) into an ordered set of bins based on its total molecule count. Genes are sorted in descending order of counts and greedily grouped so that each bin holds a roughly equal number of molecules (the most abundant genes therefore end up alone or in small bins, while many low-count genes share a bin). The resulting gene → bin assignment is written as a JSON file that can be consumed by [split-mol2bin](split_mol2bin.md) to split a molecule-level TSV into per-bin files.

This is the first of the two commands that replace the (deprecated) [split-molecule-counts](split_molecule_counts.md) command. It performs only the gene → bin assignment step; the actual splitting of molecules is done separately by [split-mol2bin](split_mol2bin.md).

## Required Options

* `--feature-tsv STR` : Input TSV file containing gene-level total counts. The gene name is expected in the first column and the count in the second column.
* `--out-json STR` : Output JSON file containing the gene → bin assignment (equivalent to the old `_bin_counts.json` output).

## Additional Options

* `--bin-count INT` : Number of bins to split the features into (roughly equal number of molecules per bin). Fewer bins may be produced if there are not enough genes. (Default: 50)
* `--in-feature-tsv-delim STR` : Delimiter for the input feature TSV file. (Default: `\t`)
* `--strip-comment-char STR` : Character that marks comment lines to be skipped at the beginning of the input file. (Default: `#`)

## Expected Output

* `[out-json]` : A JSON file containing a list of objects, one per gene, in the form:

```json
[{"gene":"GeneA","count":100,"bin":1},{"gene":"GeneB","count":80,"bin":2}, ...]
```

  The `bin` field is 1-based; a value of `0` indicates that the gene was not assigned to any bin.

## Full Usage

The full usage of `spatula assign-feature2bin` can be viewed with the `--help` option:

```
$ ./spatula assign-feature2bin --help
[./spatula assign-feature2bin] -- Assign features (genes) into ordered bins based on total counts, writing the bin assignment as JSON

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --feature-tsv          [STR: ]             : TSV file containing gene-level total counts (gene name in column 1, count in column 2)
   --out-json             [STR: ]             : Output JSON file containing the gene->bin assignment (equivalent to the old _bin_counts.json)

== Key Parameters ==
   --bin-count            [INT: 50]           : Number of bins to split the features into (roughly equal number of molecules per bin)

== Auxilary Input/Output Parameters ==
   --in-feature-tsv-delim [STR: 	]            : Delimiter for the input feature TSV file
   --strip-comment-char   [STR: #]            : Character to strip from the beginning of comment lines in the input file (if any)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
