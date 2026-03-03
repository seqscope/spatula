# spatula append_topk_tsv

## Summary 

`spatula append-topk-tsv` modifies FICTURE2 output of `punkst lda4hex` command by appending `topK` and `topP` columns with options to reorder the columns. This tool makes it easier to analyze the FICTURE2 output for certain downstream applications used in the `cartloader` pipeline.

* Input: Takes FICTURE2 output of `punkst lda4hex` command, which contains a few header columns (e.g. `random_key`, `X`, `Y`, or `cell_id`) followed by the posterior probabilities of individual topics from a LDA model. The input can also be used for other purposes, as long as each numeric part of the row represents posterior probabilities that sum to 1. Additionally, a JSON file containing the metadata can be provided.
* Output: Modified TSV file with `topK` and `topP` columns appended. Optionally, the columns can be reordered based on the total count of each topic across all cells.


An example usage of the tool is as follows:

```sh
## Usage with metadata JSON file -- add topK and topP columns and remove random_key
spatula append-topk-tsv --in-tsv /path/to/input.tsv.gz --in-json /path/to/input.json \
                    --out /path/to/output.tsv.gz --reorder
```

See below for a more detailed usage description.

## Required Options
* `--in-tsv STR`: Input TSV file containing header columns followed by topic probabilities.
* `--out-tsv STR`: Output TSV file with `topK` and `topP` columns appended with additional modifications as specified.

## Options for Reordering Columns
* `--reorder` : If turned ON, the columns in the output file are reordered based on the total count of each topic across all cells. This requires `--in-model` and `--out-model` options to be provided as well.
* `--in-model` : The model file containing gene (feature) name followed by the total counts per topic in pseudobulk format. This is required when `--reorder` is turned ON.
* `--out-model` : The output model file with the columns reordered based on the total counts per topic that will be generated as the result of `--reorder` option.

## Additional Options
* `--in-json STR`: Input JSON file containing metadata of `sptsv` format. If provided, column indices for random key and other header names are inferred from the JSON file.
* `--keep-random-key` : If turned ON, the `random_key` column is retained in the output file. By default, this column is removed.
* `--icol-random-key INT` : Column index (0-based) for the `random_key` column in the input TSV file. Default is -1, which indicates that the column is not present. If `--in-json` is provided, this option will be ignored.
* `--offset-model INT` : Column index (0-based) that indicates the beginning of the counts per topic in the input model file. Default is 1, which assumes that the first column is the gene (feature) name.
* `--offset-data INT` : Column index (0-based) that indicates the beginning of the topic probabilities in the input TSV file. Default is 3, which assumes that the first three columns are header columns (e.g. `random_key`, `X`, `Y`).
* `--topK STR` : Column name for the `topK` column to be appended. Default is `topK`.
* `--topP STR` : Column name for the `topP` column to be appended. Default is `topP`.

## Expected Output

* The output TSV file specified by `--out-tsv` will contain the same header columns as the input TSV file, followed by the topic probabilities, and then the appended `topK` and `topP` columns. If `--reorder` is turned ON, the topic columns will be reordered based on the total counts per topic.
* The output model file specified by `--out-model` will contain the same gene (feature) names as the input model file, followed by the reordered total counts per topic if `--reorder` is turned ON.

## Full Usage 

The full usage of `spatula append-topk-tsv` can be viewed with the `--help` option:

```
$ ./spatula append-topk-tsv --help
[./spatula append-topk-tsv] -- Append topK and topP columns to the input TSV file

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --in-tsv          [STR: ]             : Input pseudobulk TSV file
   --in-json         [STR: ]             : Input JSON file
   --in-model        [STR: ]             : Input model file (required with --reorder)
   --out-tsv         [STR: ]             : Output TSV file
   --out-model       [STR: ]             : Output model file (required with --reorder)
   --reorder         [FLG: OFF]          : Reorder the columns in the output file based on the total count
   --keep-random-key [FLG: OFF]          : Keep the random key in the output file

== Expected columns in input and output ==
   --icol-random-key [INT: -1]           : Column index for the random key, -1 if not present
   --offset-model    [INT: 1]            : Column index for the beginning of the input model file
   --offset-data     [INT: 3]            : Column index for the beginning of the input tsv file
   --topK            [STR: topK]         : Column name for topK
   --topP            [STR: topP]         : Column name for topP


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```