# spatula draw-tsv

## Summary

`spatula draw-tsv` draws an image of TSV-formatted spatial gene expression data, coloring each transcript according to a per-gene color map. Given a generic transcript TSV (with X/Y/gene/count columns) and a color map file that assigns an RGB hex color to each gene, it accumulates colored intensity per pixel to produce an RGB image of the spatial distribution of expression.

This is the TSV-input counterpart of [draw-sge](draw_sge.md) (which reads an SGE matrix).

!!! note
    A `--minmax` bounding box (with `xmin`/`xmax`/`ymin`/`ymax`) defines the drawn extent.

## Required Options

* `--tsv STR` : Input TSV file.
* `--cmap STR` : Color map file containing gene name and RGB color.
* `--out STR` : Output file name.

## Additional Options

### Input
* `--minmax STR` : Bounding box information. Expects `xmin`/`xmax`/`ymin`/`ymax`.

### Column Names
* `--tsv-colname-x STR` : Column name for the X coordinate in the input TSV. (Default: `X`)
* `--tsv-colname-y STR` : Column name for the Y coordinate in the input TSV. (Default: `Y`)
* `--tsv-colname-gene STR` : Column name for the gene in the input TSV. (Default: `gene`)
* `--tsv-colname-count STR` : Column name for the count in the input TSV. (Default: `count`)
* `--cmap-colname-gene STR` : Column name for the gene in the color map file. (Default: `gene`)
* `--cmap-colname-hex STR` : Column name for the hex color in the color map file. (Default: `hex`)
* `--skip-tsv-header` : Skip the header in the input TSV file.

### Output
* `--coord-per-pixel FLT` : Number of coordinate units per pixel. (Default: `1.00`)
* `--intensity-per-count FLT` : Intensity per count in RGB color space. (Default: `0.10`)

## Expected Output

* `[out]` : A PNG image visualizing the colored spatial expression.

## Full Usage

The full usage of `spatula draw-tsv` can be viewed with the `--help` option:

```
$ ./spatula draw-tsv --help
[./spatula draw-tsv] -- Draw the image of TSV-formatted spatial gene expression data

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input files ==
   --minmax              [STR: ]             : Bounding box information. Expects xmin/xmax/ymin/ymax
   --tsv                 [STR: ]             : Input TSV file
   --cmap                [STR: ]             : Color map file containing gene name and RGB color

== Column names ==
   --tsv-colname-x       [STR: X]            : Column name for the X coordinate in the input TSV file
   --tsv-colname-y       [STR: Y]            : Column name for the Y coordinate in the input TSV file
   --tsv-colname-gene    [STR: gene]         : Column name for the gene in the input TSV file
   --tsv-colname-count   [STR: count]        : Column name for the count in the input TSV file
   --cmap-colname-gene   [STR: gene]         : Column name for the gene in the color map file
   --cmap-colname-hex    [STR: hex]          : Column name for the hex color in the color map file
   --skip-tsv-header     [FLG: OFF]          : Skip the header in the input TSV file

== Output options ==
   --coord-per-pixel     [FLT: 1.00]         : Number of coordinate units per pixel
   --intensity-per-count [FLT: 0.10]         : Intensity per count in RGB color space

== Output Options ==
   --out                 [STR: ]             : Output file name


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
