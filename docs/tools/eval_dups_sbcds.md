# spatula eval-dups-sbcds

## Summary 

`spatula eval-dups-sbcds` examines duplicates in the spatial barcode dictionary
generated by `spatula build-sbcds` and produces summary statistics.

Here is a summary of the main features:

* Input: Takes a spatial barcode directory from `build-sbcds`.
* Output: In the output directory, it produces (a) a complete list of duplicate barcodes, with the number of occurrences and the number of unique tiles observed, (b) per-tile summary metrics of the number of duplicate barcodes, and (c) the histogram of the observed counts of each duplicate.

An example use of the tool is as follows:

```sh
spatula eval-dups-sbcds --sbcd /path/to/sbcd/dir/ \
                        --out  /path/to/output_prefix 
```
See below for a more detailed usage description.

## Required options

* `--sbcd` : The path to the directory containing spatial barcode files, generated by `build-sbcds`.
* `--out` : The prefix for the output files. See [Expected Output](#expected-output) for more details.

## Additional options

* `--match-len` : The length of the spatial barcode to consider for a match. The default  is 27, and the maximum possible value is 27.

## Expected Output

With `[out_prefix]` as the prefix, the following files will be created:

- `[out_prefix].dups.sbcds.tsv.gz` : A compressed tsv file containing the following entries in each line:
    - `barcode` : Spatial barcode that is observed multiple times across the tiles.
    - `ntiles` : Number of unique tiles where the spatial barcodes were observed.
    - `dupcount` : The number of occurrences of the spatial barcode across the tiles.
- `[out_prefix].tiles.tsv` : A tab-delimited file summarizing the duplicate barcode statistics. It contains the following entries in each row:
    - `tile` : The ID of the tile, in the format of `lane_tile`.
    - `total` : The total number of spatial barcodes found in the tile.
    - `uniq` : The number of unique spatial barcodes (across tiles) found in the tile.
    - `dups_uniq` : The number of non-unique spatial barcodes in the tile, counting only once per tile.
    - `dups_within` : The number of non-unique spatial barcodes observed only within a single tile, counting only once per tile.
    - `dups_between` : The number of non-unique spatial barcodes observed across multiple tiles, counting only once per tile.
- `[out_prefix].hist.tsv` : The histogram of the observed counts of each duplicate, across all tiles. It contains the following entries in each row.
    - `dupcount` : The number of occurrences of the spatial barcode across the tiles.
    - `num` : The number of spatial barcodes observed with the `dupcount` occurrences.


## Full Usage 

The full usage of the software tool is as follows:

```sh
$ ./spatula eval-dups-sbcds --help
[./spatula eval-dups-sbcds] -- Evaluate duplicates in spatial barcodes

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --sbcd      [STR: ]             : Spatial barcode dictionary generated from 'build-sbcds' command
   --match-len [INT: 27]           : Length of HDMI spatial barcodes to require perfect matches

== Output Options ==
   --out       [STR: ]             : Prefix of output files (.dups.sbcds.tsv.gz, .tiles.tsv, .hist.tsv)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```