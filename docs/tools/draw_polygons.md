# spatula draw-polygons

## Summary

`spatula draw-polygons` draws polygons from a GeoJSON input, rendering an image in which each polygon is filled with a color determined by one of its properties. The color can come from a cluster attribute mapped through an RGB color-map file, or the cluster names can themselves be interpreted directly as `RRGGBB` hex values. This is useful for visualizing cell/segment boundaries or cluster assignments produced elsewhere in the pipeline.

See [pixels2cells](pixels2cells.md) and [tsv-add-cell-id](tsv_add_cell_id.md) for producing polygon/cell assignments.

## Required Options

* `--boundary STR` : GeoJSON file to draw polygons from.
* `--out STR` : Output file name.

## Additional Options

### Input
* `--ullr STR` : Upper-left and lower-right coordinates as `ulx,uly,lrx,lry`.
* `--range STR` : TSV file specifying the boundary box.
* `--clust-attr STR` : Cluster (property) name used to determine the color mapping. (Default: `topK`)
* `--cmap STR` : RGB color mapping file.
* `--cluster-as-rgb` : Interpret cluster names directly as RGB values in `RRGGBB` format.

### Output
* `--scale FLT` : Scale factor — pixels per coordinate unit. (Default: `1.00`)
* `--opacity FLT` : Opacity of the polygons, from `0.0` to `1.0`. (Default: `0.50`)

## Expected Output

* `[out]` : An image with the rendered polygons.

## Full Usage

The full usage of `spatula draw-polygons` can be viewed with the `--help` option:

```
$ ./spatula draw-polygons --help
[./spatula draw-polygons] -- Draw polygons based on JSON input

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --boundary       [STR: ]             : GeoJSON file to draw polygons from
   --ullr           [STR: ]             : Specify upper-left and lower-right coordinates as 'ulx,uly,lrx,lry'
   --range          [STR: ]             : TSV file specifying the boundary box
   --clust-attr     [STR: topK]         : Cluster name to determine color mapping
   --cmap           [STR: ]             : RGB color mapping file
   --cluster-as-rgb [FLG: OFF]          : Indicate that cluster names are RGB values as 'RRGGBB'

== Output Options ==
   --scale          [FLT: 1.00]         : Scale factor: pixels per coordinate unit
   --out            [STR: ]             : Output file name
   --opacity        [FLT: 0.50]         : Opacity of the polygons (0.0 to 1.0)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
