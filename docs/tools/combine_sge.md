# spatula combine-sge

## Summary

`spatula combine-sge` combines multiple spatial gene expression (SGE) matrices into a single SGE, arranging the individual SGEs in space according to a layout file. Each input SGE contributes its barcodes, features, and matrix, and the tool re-tiles them either by row/column grid position or by explicit x/y offsets. This is useful for stitching together per-tile or per-section SGE outputs into one dataset.

See [dge2sge](dge2sge.md) for producing individual SGE matrices, and [convert-sge](convert_sge.md) for converting an SGE to a generic TSV.

## Required Options

* `--layout STR` : Layout file, with one row per input SGE directory. Each row contains `[sgedir]` and either `[row]`/`[col]` (grid mode) or `[x_offset]`/`[y_offset]` columns.
* `--out STR` : Output directory.

## Additional Options

* `--bcd STR` : Shared barcode file name within each SGE directory. (Default: `barcodes.tsv.gz`)
* `--ftr STR` : Shared feature file name within each SGE directory. (Default: `features.tsv.gz`)
* `--mtx STR` : Shared matrix file name within each SGE directory. (Default: `matrix.mtx.gz`)
* `--minmax STR` : Shared `minmax.tsv` file name within each SGE directory. Required in row/col mode. (Default: `barcodes.minmax.tsv`)
* `--out-minmax-fixed` : Do not update the output minmax coordinates based on the observed points.
* `--mode-rowcol` : Activate row/col (grid) mode instead of offset mode. (Default: OFF)

## Expected Output

* `[out]/` : Output directory containing the combined SGE (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`, and associated metadata).

## Full Usage

The full usage of `spatula combine-sge` can be viewed with the `--help` option:

```
$ ./spatula combine-sge --help
[./spatula combine-sge] -- Combine multiple SGE files

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --layout           [STR: ]             : Layout file, each containing [sgedir] and [row]/[col] or [x_offset]/[y_offset] as columns
   --bcd              [STR: barcodes.tsv.gz] : Shared barcode file path (e.g. barcodes.tsv.gz)
   --ftr              [STR: features.tsv.gz] : Shared feature file path (e.g. feature.tsv.gz)
   --mtx              [STR: matrix.mtx.gz] : Shared matrix file path (e.g. matrix.mtx.gz)
   --minmax           [STR: barcodes.minmax.tsv] : Shared minmax.tsv file path (e.g. barcodes.minmax.tsv) - required in [row]/[col] mode
   --out-minmax-fixed [FLG: OFF]          : Do not update output minmax coordinates based on the observed points
   --mode-rowcol      [FLG: OFF]          : Activate rowcol mode (default: false)

== Output Options ==
   --out              [STR: ]             : Output directory


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
