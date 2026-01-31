# spatula join-pixel-decode

## Summary 

`spatula join-pixel-decode` matches individual molecules (from a "molecular" TSV file) to pixel-level decoding results (from "decoded" TSV files, e.g., output from FICTURE) based on spatial proximity. For each molecule, it searches for the nearest pixel in the decoded file(s) within a specified distance and augments the molecular data with the pixel's attributes (e.g., cluster assignment, posterior probabilities).

## Required Options

* `--mol-tsv STR` : Input TSV file containing individual molecules with X and Y coordinates.
* `--decode-prefix-tsv STR` : Argument specifying a pixel-level output file in the format `[prefix],[tsv-path]`. This option can be used multiple times to join multiple decoded files.
* `--out-prefix STR` : Output prefix for the joined TSV file.

## Additional Options

* `--colname-mol-x STR` : Column name for the X coordinate in the molecular TSV. (Default: X)
* `--colname-mol-y STR` : Column name for the Y coordinate in the molecular TSV. (Default: Y)
* `--colname-decode-x STR` : Column name for the X coordinate in the decoded TSV. (Default: x)
* `--colname-decode-y STR` : Column name for the Y coordinate in the decoded TSV. (Default: y)
* `--tile-size FLT` : Size of the tiles (square regions) used for spatial indexing and parallel processing, in the same units as coordinates. (Default: 500.00)
* `--max-dist FLT` : Maximum distance allowed to link a molecule to a pixel. If the nearest pixel is further than this distance, values are set to NA. (Default: 0.50)
* `--precision INT` : Precision (number of decimal places) for writing coordinate values in the output. (Default: 3)
* `--mu-scale FLT` : Scaling factor for the coordinate resolution. (Default: 1.00)
* `--threads INT` : Number of threads to use for parallel processing. (Default: 1)
* `--colnames-include STR` : Comma-separated list of columns from the molecular TSV to include in the output.
* `--colnames-exclude STR` : Comma-separated list of columns from the molecular TSV to exclude from the output.
* `--out-max-k INT` : Maximum number of pixel-level factors (e.g., top cluster indices) to include from the decoded file. (Default: 1)
* `--out-max-p INT` : Maximum number of posterior probabilities to include from the decoded file. (Default: 1)
* `--tmp-dir STR` : Directory to store temporary intermediate files.
* `--out-suffix-tsv STR` : Suffix for the output filename. (Default: .tsv)

## Expected Output

The tool produces a single TSV file named `[out-prefix].tsv` (or specified suffix). The file contains:
1.  **Coordinates**: `X` and `Y` columns from the molecular file.
2.  **Molecular Attributes**: Columns retained from the input molecular TSV (based on include/exclude options).
3.  **Joined Attributes**: For each input `decode-prefix-tsv` (with prefix `P`), columns such as `P[ColName]` are added. Common columns include `PK1` (top cluster), `PP1` (max probability), etc., corresponding to the nearest pixel in that decoded file. If no pixel is found within `--max-dist`, these columns will contain `NA`.

## Full Usage 

The full usage of `spatula join-pixel-decode` can be viewed with the `--help` option:

```
$ ./spatula join-pixel-decode --help   
[./spatula join-pixel-decode] -- Join pixel-level-decode output from FICTURE2 with raw transcript-level TSV files

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --mol-tsv           [STR: ]             : TSV file containing individual molecules
   --decode-prefix-tsv [V_STR: ]           : TSV file containing pixel-level factors in [prefix],[tsv-path] format
   --out-prefix        [STR: ]             : Output prefix for the joined TSV files
   --tmp-dir           [STR: ]             : Temporary directory for intermediate files

== Key Parameters ==
   --tile-size         [FLT: 500.00]       : Tile size to create temporary files for binning
   --max-dist          [FLT: 0.50]         : Maximum distance in um to consider a match
   --mu-scale          [FLT: 1.00]         : Scale factor for the resolution - divide by mu_scale in the output
   --precision         [FLT: 2.1e-314]     : Output precision below the decimal point
   --threads           [INT: 1]            : Number of threads to use for processing

== Expected columns in input and output ==
   --colname-mol-x     [STR: X]            : Column name for X-axis for molecular TSV
   --colname-mol-y     [STR: Y]            : Column name for Y-axis for molecular TSV
   --colname-decode-x  [STR: X]            : Column name for X-axis for decoded TSV
   --colname-decode-y  [STR: Y]            : Column name for Y-axis for decoded TSV
   --colnames-include  [STR: ]             : Comma-separated column names to include in the output TSV file
   --colnames-exclude  [STR: ]             : Comma-separated column names to exclude in the output TSV file
   --out-max-k         [INT: 1]            : Maximum number of pixel-level factors to include in the joined output. (Default : 1)
   --out-max-p         [INT: 1]            : Maximum number of pixel-level posterior probabilities to include in the joined output. (Default : 1)

== Output File suffixes ==
   --out-suffix-tsv    [STR: .tsv]         : Suffix for the output TSV file


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```