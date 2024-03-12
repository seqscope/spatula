# spatula draw-3way

## Summary 

`spatula draw-3way` is a tool to visualize the spatial distribution of 1st-seq and 2nd-seq reads of Seq-Scope data for each 'Chip'.

**IMPORTANT** [ImageMagick](https://imagemagick.org/script/download.php) must be installed to use this tool.

A typical use case is as follows:

* Input: Takes (a) an output `combine-sbcds`, (b) (multiple) outputs `match-sbcds`, and (c) output from `dge2sge` commands from individual tools.
* Output: Produces an 2D image that plots the (a) all 1st-seq spatial barcodes from `combine-sbcds` in blue, (b) all 2nd-seq spatial barcodes from `match-sbcds` in green, and (c) all 2nd-seq spatial barcodes aligned to genes from `dge2sge` in red.

A typical example is as follows:

```sh
spatula draw-3way --manifest /path/to/combine/sbcds/output/dir/manifest.tsv \
                  --nbcd     /path/to/combine/sbcds/output/dir/1_1.sbcds.sorted.tsv.gz \
                  --nmatch   /path/to/match/sbcds/output/prefix1.match.sorted.uniq.tsv.gz \
                  --nmatch   /path/to/match/sbcds/output/prefix2.match.sorted.uniq.tsv.gz \
                  --ngebcd   /path/to/dge2sge/output/dir/barcodes.tsv.gz \
                  --out      /path/to/output/image.png                
```


See below for a more detailed usage description.

## Required options
* `--manifest` : The `manifest.tsv` file from the `combine-sbcds` file that contains the summary of the spatial coordinate of a Seq-Scope Chip.
* `--nbcd` : The spatial barcode file from the `combine-sbcds` command. The filename is usually `1_1.sbcds.sorted.tsv.gz`.
* `--nmatch` : (Multiple allowed) The spatial barcode file from the `match-sbcds` command. The filename is usually `prefix.match.sorted.uniq.tsv.gz`.
* `--ngebcd` : The spatial barcode file from the `dge2sge` command. The filename is usually `barcodes.tsv.gz`.
* `--out` : The output filename of the image. Currently, `.png` is supported. 

## Additional Options

* `--coord-per-pixel` : The number of coordinates to be collapsed into a pixel as a factor to divide the input coordinate with. The default is 1000.0.
* `--color-nbcd` : The RGB hex color code for the 1st-seq spatial barcodes in the order of G,B,R. The default is `#004800` (blue).
* `--color-nmatch` : The RGB hex color code for the 2nd-seq spatial barcodes in the order of G,B,R. The default is `#640000` (green).
* `--color-nge` : The RGB hex color code for the 2nd-seq spatial barcodes aligned to genes in the order of G,B,R. The default is `#000032` (red).
* `--icol-x-nbcd` : The (0-based) column index of the X coordinate in the 1st-seq spatial barcode file. The default is 3.
* `--icol-y-nbcd` : The (0-based) column index of the Y coordinate in the 1st-seq spatial barcode file. The default is 4.
* `--icol-x-nmatch` : The (0-based) column index of the X coordinate in the 2nd-seq spatial barcode file. The default is 3.
* `--icol-y-nmatch` : The (0-based) column index of the Y coordinate in the 2nd-seq spatial barcode file. The default is 4.
* `--icol-x-ngebcd` : The (0-based) column index of the X coordinate in the 2nd-seq spatial barcode file. The default is 5.
* `--icol-y-ngebcd` : The (0-based) column index of the Y coordinate in the 2nd-seq spatial barcode file. The default is 6.
* `--icol-gene-ngebcd` : The (0-based) column index of the gene count in the 2nd-seq spatial barcode file. The default is 7.
* `--isubcol-gene-ngebcd` : The (0-based) sub column index of the gene count (among the comma-separated-values) in the 2nd-seq spatial barcode file. The default is 0.
* `--id-manifest` : The ID of the tile in the manifest file to use. The default is `1_1`.

## Expected Output

The output `[out]` will be created as a PNG file containing the image of the input points.

## Full Usage 

The full usage of `spatula draw-3way` can be viewed with the `--help` option:

```
$ ./spatula draw-3way --help
[./spatula draw-3way] -- Draw the 3-way image from the output of sttools pipeline

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input files ==
   --manifest            [STR: ]             : Manifest file from combine-sbcds. xmin/xmax/ymin/ymax will be automatically detected
   --nbcd                [STR: ]             : Spatial barcode dictionary generated from 'combine-sbcds' command
   --nmatch              [V_STR: ]           : Spatial barcode dictionary generated from 'match-sbcds' command
   --ngebcd              [STR: ]             : Spatial barcode dictionary generated from alignment pipeline

== Input options ==
   --icol-x-nbcd         [INT: 3]            : 0-based index of the column for x in nbcd
   --icol-y-nbcd         [INT: 4]            : 0-based index of the column for y in nbcd
   --icol-x-nmatch       [INT: 3]            : 0-based index of the column for x in nmatch
   --icol-y-nmatch       [INT: 4]            : 0-based index of the column for y in nmatch
   --icol-x-ngebcd       [INT: 5]            : 0-based index of the column for x in ngebcd
   --icol-y-ngebcd       [INT: 6]            : 0-based index of the column for y in ngebcd
   --icol-gene-ngebcd    [INT: 7]            : 0-based index of the column for gene count in ngebcd
   --isubcol-gene-ngebcd [INT: 0]            : 0-based index of the sub column for gene count in ngebcd (0: Gene, 1: GeneFull, ...)
   --id-manifest         [STR: 1_1]          : ID of the tile in the manifest file to use

== Output options ==
   --coord-per-pixel     [FLT: 1000.00]      : Number of coordinate units per pixel
   --color-nbcd          [STR: #004800]      : RGB hex color code for nbcd per observation
   --color-nmatch        [STR: #640000]      : RGB hex color code for nmatch per observation
   --color-nge           [STR: #000032]      : RGB hex color code for nge per observation

== Output Options ==
   --out                 [STR: ]             : Output file name


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```