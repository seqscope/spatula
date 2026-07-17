# spatula filter-tsv

## Summary

`spatula filter-tsv` filters a generic spatial transcript TSV file by spatial coordinates and/or by gene identity. Records can be restricted to a bounding box (`--xmin`/`--xmax`/`--ymin`/`--ymax`) or to the interior of a polygon (`--polygon`, GeoJSON), and/or limited to a list of gene names (`--filt-gene`) or gene IDs (`--filt-gid`). Alongside the filtered transcript TSV, it writes an updated feature count file and a minmax coordinate file.

## Required Options

* `--in-tsv STR` : Input (unsorted) TSV file.
* `--out-prefix STR` : Prefix of output files.

## Additional Options

### Output Suffixes
* `--tsv-suffix STR` : Suffix for the output transcript TSV file. (Default: `.transcripts.tsv.gz`)
* `--ftr-suffix STR` : Suffix for the output feature file. (Default: `.features.tsv.gz`)
* `--minmax-suffix STR` : Suffix for the output minmax file. (Default: `.features.tsv.gz`)

### Column Names
* `--colname-x STR` : Column name for the X coordinate. (Default: `X`)
* `--colname-y STR` : Column name for the Y coordinate. (Default: `Y`)
* `--colname-gene STR` : Column name for the gene name. (Default: `gene`)
* `--colname-gid STR` : Column name for the gene ID. (Default: `gene_id`)
* `--colname-count STR` : Column name for the count. (Default: `count`)

### Spatial Filtering
* `--xmin FLT` / `--xmax FLT` : Minimum/maximum X value. (Defaults: `-inf` / `inf`)
* `--ymin FLT` / `--ymax FLT` : Minimum/maximum Y value. (Defaults: `-inf` / `inf`)
* `--polygon STR` : GeoJSON file for polygon-based filtering.

### Gene Filtering
* `--filt-gene STR` : Include only gene names present in the given list file.
* `--filt-gid STR` : Include only gene IDs present in the given list file.

## Expected Output

* `[out-prefix][tsv-suffix]` : The filtered transcript TSV file.
* `[out-prefix][ftr-suffix]` : Updated feature counts for the filtered data.
* `[out-prefix][minmax-suffix]` : Minmax coordinate file for the filtered data.

## Full Usage

The full usage of `spatula filter-tsv` can be viewed with the `--help` option:

```
$ ./spatula filter-tsv --help
[./spatula filter-tsv] -- Filter TSV based on spatial coordinates or gene list

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --in-tsv        [STR: ]             : Input (unsorted) TSV file
   --out-prefix    [STR: ]             : Prefix of output files
   --tsv-suffix    [STR: .transcripts.tsv.gz] : Suffix for output TSV file
   --ftr-suffix    [STR: .features.tsv.gz] : Suffix for output feature file
   --minmax-suffix [STR: .features.tsv.gz] : Suffix for output minmax file

== Column Names ==
   --colname-x     [STR: X]            : Column name for X-axis
   --colname-y     [STR: Y]            : Column name for Y-axis
   --colname-gene  [STR: gene]         : Column name for gene name
   --colname-gid   [STR: gene_id]      : Column name for gene ID
   --colname-count [STR: count]        : Column name for count

== Spatial Filtering options ==
   --xmin          [FLT: -inf]         : Minimum x-axis value
   --xmax          [FLT: inf]          : Maximum x-axis value
   --ymin          [FLT: -inf]         : Minimum y-axis value
   --ymax          [FLT: inf]          : Maximum y-axis value
   --polygon       [STR: ]             : GeoJSON file for polygon-based filtering

== Gene Filtering options ==
   --filt-gene     [STR: ]             : Only Include gene names present in the list
   --filt-gid      [STR: ]             : Only Include gene ids present in the list


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
