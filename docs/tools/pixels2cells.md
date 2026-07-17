# spatula pixels2cells

## Summary

`spatula pixels2cells` assigns pixels (individual molecule positions) to cell IDs based on cell boundary polygons provided as a GeoJSON file. Each pixel's `X`/`Y` coordinate is tested against the polygons, and the matching cell ID is written to a new column. Pixels that fall inside no polygon are labeled with the unassigned cell ID. The assignment uses a rasterized grid whose resolution is controlled by `--um-per-grid`.

For an exact, multithreaded point-in-polygon alternative that also supports boundary expansion, see [tsv-add-cell-id](tsv_add_cell_id.md).

## Required Options

* `--boundary STR` : GeoJSON file containing cell IDs and polygons.
* `--pixel STR` : Pixel-level TSV data containing individual molecule positions.
* `--out STR` : Output TSV file name with cell assignments.

## Additional Options

### GeoJSON Input
* `--cell-id-attr STR` : Attribute name for cell IDs in the GeoJSON properties. (Default: `cell_id`)

### Pixel TSV Input
* `--colname-x STR` : Column name for the pixel X coordinate. (Default: `X`)
* `--colname-y STR` : Column name for the pixel Y coordinate. (Default: `Y`)

### Output
* `--out-colname-cell-id STR` : Column name for the cell ID in the output TSV. (Default: `cell_id`)
* `--unassigned-cell-id STR` : Cell ID assigned to pixels outside all polygons. (Default: `UNASSIGNED`)

### Grid / Bounding Box
* `--um-per-grid FLT` : Grid size in microns. (Default: `1.0`)
* `--range STR` : TSV file specifying `xmin`, `ymin`, `xmax`, `ymax` for the area of interest.
* `--ullr STR` : Upper-left and lower-right coordinates as `ulx,uly,lrx,lry`.

!!! note
    If both `--range` and `--ullr` are omitted, the boundary file is read twice in order to determine the bounding box.

## Expected Output

* `[out]` : The input pixel TSV with an added cell-ID column.

## Full Usage

The full usage of `spatula pixels2cells` can be viewed with the `--help` option:

```
$ ./spatula pixels2cells --help
[./spatula pixels2cells] -- Assign pixels to cell IDs based on boundary polygons

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input GeoJSON options ==
   --boundary            [STR: ]             : GeoJSON file containing cell IDs and polygons
   --cell-id-attr        [STR: cell_id]      : Attribute name for cell IDs in the GeoJSON properties

== Input Pixel TSV options ==
   --pixel               [STR: ]             : Pixel-level TSV data containing individual molecule positions
   --colname-x           [STR: X]            : Column name for pixel X coordinate
   --colname-y           [STR: Y]            : Column name for pixel Y coordinate

== Output Options ==
   --out                 [STR: ]             : Output TSV file name with cell assignments
   --out-colname-cell-id [STR: cell_id]      : Column name for cell ID in the output TSV
   --unassigned-cell-id  [STR: UNASSIGNED]   : Cell ID to assign for unassigned pixels

== Other Input Options ==
   --um-per-grid         [FLT: 1.00]         : Grid size in microns (default: 1.0)
   --range               [STR: ]             : TSV file specifying xmin, ymin, xmax, ymax for the area of interest. If both --range and --ullr is absent, the boundary file will be read twice to determine the bounding box.
   --ullr                [STR: ]             : Specify upper-left and lower-right coordinates as 'ulx,uly,lrx,lry'. If both --range and --ullr is absent, the boundary file will be read twice to determine the bounding box.


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
