# spatula paste-pixel-tsv

## Summary

`spatula paste-pixel-tsv` pastes pixel-level factor output from FICTURE2 together with raw transcript-level TSV files, matching records by their `X`/`Y` coordinates. It is used to append the pixel-level factor assignments and posterior probabilities to the raw transcript TSV so that downstream tools can operate on a single combined table. Multiple pixel factor files can be supplied and combined.

This tool is the FICTURE2 counterpart of [join-pixel-tsv](join_pixel_tsv.md); see also [join-pixel-decode](join_pixel_decode.md).

## Required Options

* `--pix-prefix-tsv V_STR` : TSV file(s) containing pixel-level factors. May be specified multiple times.
* `--out-tsv STR` : Output TSV file.

## Additional Options

* `--colname-x STR` : Column name for the X coordinate. (Default: `X`)
* `--colname-y STR` : Column name for the Y coordinate. (Default: `Y`)
* `--colnames-include STR` : Comma-separated list of column names to include in the output TSV.
* `--colnames-exclude STR` : Comma-separated list of column names to exclude from the output TSV.
* `--out-max-k INT` : Maximum number of pixel-level factors to include in the joined output. (Default: 1)
* `--out-max-p INT` : Maximum number of pixel-level posterior probabilities to include in the joined output. (Default: 1)

## Expected Output

* `[out-tsv]` : The combined TSV file with pixel-level factor and posterior-probability columns pasted onto the transcript records.

## Full Usage

The full usage of `spatula paste-pixel-tsv` can be viewed with the `--help` option:

```
$ ./spatula paste-pixel-tsv --help
[./spatula paste-pixel-tsv] -- Paste pixel-level output from FICTURE2 with raw transcript-level TSV files

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --pix-prefix-tsv   [V_STR: ]           : TSV file containing pixel-level factors
   --out-tsv          [STR: ]             : Output TSV file

== Expected columns in input and output ==
   --colname-x        [STR: X]            : Column name for X-axis
   --colname-y        [STR: Y]            : Column name for Y-axis
   --colnames-include [STR: ]             : Comma-separated column names to include in the output TSV file
   --colnames-exclude [STR: ]             : Comma-separated column names to exclude in the output TSV file
   --out-max-k        [INT: 1]            : Maximum number of pixel-level factors to include in the joined output. (Default : 1)
   --out-max-p        [INT: 1]            : Maximum number of pixel-level posterior probabilities to include in the joined output. (Default : 1)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
