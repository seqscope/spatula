# spatula tsv-add-cell-id

## Summary

`spatula tsv-add-cell-id` annotates a transcript-level TSV with cell IDs derived from cell boundary polygons, using exact point-in-polygon tests and multithreaded processing. Boundaries can be supplied either as a GeoJSON file or as a CSV of polygon vertices. Each transcript's `X`/`Y` coordinate is tested against the polygons and the matching cell ID is appended as a new column; transcripts outside all polygons receive the unassigned value. Boundaries can optionally be expanded outward so that nearby transcripts are assigned to the closest cell.

This is a faster, exact alternative to the grid-based [pixels2cells](pixels2cells.md).

## Required Options

* `--tsv STR` : Transcript-level TSV(.gz) file with `X`/`Y` coordinates (streamed).
* `--out STR` : Output TSV(.gz) file with the added cell-ID column.
* One boundary source is required — either `--boundaries-geojson` **or** `--boundaries-csv` (see below).

## Additional Options

### Input Transcript TSV
* `--colname-x STR` : Column name for the X coordinate in the TSV. (Default: `X`)
* `--colname-y STR` : Column name for the Y coordinate in the TSV. (Default: `Y`)

### Input Boundaries (one required)
* `--boundaries-geojson STR` : GeoJSON file with cell boundary polygons.
* `--boundaries-csv STR` : CSV(.gz) with `cell_id`, `vertex_x`, `vertex_y` columns.
* `--cell-id-attr STR` : GeoJSON property name holding the cell ID. (Default: `cell_id`)
* `--csv-colname-cell-id STR` : CSV column name for cell ID. (Default: `cell_id`)
* `--csv-colname-x STR` : CSV column name for vertex X. (Default: `vertex_x`)
* `--csv-colname-y STR` : CSV column name for vertex Y. (Default: `vertex_y`)

### Output
* `--out-colname-cell-id STR` : Column name for the added cell ID. (Default: `cell_id`)
* `--unassigned-cell-id STR` : Value for transcripts inside no polygon. (Default: `UNASSIGNED`)

### Boundary Expansion
* `--expand-um FLT` : Expand cell boundaries outward by this distance; a point outside all polygons is assigned to the nearest polygon within this distance (`0` = off). (Default: `0.00`)

### Performance
* `--threads INT` : Number of threads for point-in-polygon assignment. (Default: 12)
* `--batch-size INT` : Number of TSV lines processed per batch. (Default: 500000)
* `--bucket-um FLT` : Grid cell size for the spatial index (`0` = auto). (Default: `0.00`)

## Expected Output

* `[out]` : The input transcript TSV with an added cell-ID column.

## Full Usage

The full usage of `spatula tsv-add-cell-id` can be viewed with the `--help` option:

```
$ ./spatula tsv-add-cell-id --help
[./spatula tsv-add-cell-id] -- Annotate transcript TSV with cell IDs from boundary polygons (exact, multithreaded)

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Transcript TSV Options ==
   --tsv                 [STR: ]             : Transcript-level TSV(.gz) file with X/Y coordinates (streamed)
   --colname-x           [STR: X]            : Column name for X coordinate in the TSV
   --colname-y           [STR: Y]            : Column name for Y coordinate in the TSV

== Input Boundary Options (one required) ==
   --boundaries-geojson  [STR: ]             : GeoJSON file with cell boundary polygons
   --boundaries-csv      [STR: ]             : CSV(.gz) with cell_id, vertex_x, vertex_y columns
   --cell-id-attr        [STR: cell_id]      : GeoJSON property name holding the cell ID
   --csv-colname-cell-id [STR: cell_id]      : CSV column name for cell ID
   --csv-colname-x       [STR: vertex_x]     : CSV column name for vertex X
   --csv-colname-y       [STR: vertex_y]     : CSV column name for vertex Y

== Output Options ==
   --out                 [STR: ]             : Output TSV(.gz) file with the added cell-id column
   --out-colname-cell-id [STR: cell_id]      : Column name for the added cell ID
   --unassigned-cell-id  [STR: UNASSIGNED]   : Value for transcripts inside no polygon

== Boundary Expansion Options ==
   --expand-um           [FLT: 0.00]         : Expand cell boundaries outward by this distance; a point outside all polygons is assigned to the nearest polygon within this distance (0 = off)

== Performance Options ==
   --threads             [INT: 12]           : Number of threads for point-in-polygon assignment
   --batch-size          [INT: 500000]       : Number of TSV lines processed per batch
   --bucket-um           [FLT: 0.00]         : Grid cell size for the spatial index (0 = auto)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
