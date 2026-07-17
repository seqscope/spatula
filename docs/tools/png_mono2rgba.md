# spatula png-mono2rgba

## Summary

`spatula png-mono2rgba` converts a monochrome (grayscale) PNG image into an RGBA image, mapping the grayscale intensity to a chosen RGB color and computing a per-pixel alpha (transparency) channel. The alpha can be constant, linear, or quadratic in the source intensity, and pixels below/above configurable thresholds can be made fully transparent. This is typically used to turn a grayscale density/mask image into a colored, transparency-aware overlay.

## Required Options

* `--in STR` : Input PNG file (grayscale).
* `--out STR` : Output PNG file.

## Additional Options

### Color / Transparency Settings
* `--rgb STR` : RGB color applied to the monochrome image, in `RRGGBB` format. (Default: `FFFFFF`)
* `--transparent-below INT` : Pixels at or below this intensity are made transparent. (Default: 0)
* `--transparent-above INT` : Pixels at or above this intensity are made transparent. (Default: 255)
* `--invert-transparency` : Invert the transparency logic. (Default: OFF)

### Alpha Calculation
* `--linear-alpha` : Use a linear alpha calculation (instead of the default constant alpha). (Default: OFF)
* `--quadratic-alpha` : Use a quadratic alpha calculation. (Default: OFF)
* `--alpha-constant-max INT` : Maximum value for constant alpha. (Default: 255)
* `--alpha-linear-coef FLT` : Coefficient for the linear alpha calculation. (Default: `1.0`)
* `--alpha-quadratic-coef FLT` : Coefficient for the quadratic alpha calculation. (Default: `16.0`)

## Expected Output

* `[out]` : The RGBA PNG image.

## Full Usage

The full usage of `spatula png-mono2rgba` can be viewed with the `--help` option:

```
$ ./spatula png-mono2rgba --help
[./spatula png-mono2rgba] -- Convert a monochrome PNG image to RGBA format

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --in                   [STR: ]             : Input PNG file (grayscale)
   --out                  [STR: ]             : Output PNG file

== Color settings ==
   --transparent-below    [INT: 0]            : Threshold for transparent pixels (default: 0)
   --transparent-above    [INT: 255]          : Threshold for transparent pixels (default: 255)
   --rgb                  [STR: FFFFFF]       : RGB color for the monochrome image (default: #FFFFFF)
   --linear-alpha         [FLG: OFF]          : Use linear alpha calculation (default: constant)
   --quadratic-alpha      [FLG: OFF]          : Use quadratic alpha calculation (default: constant)
   --invert-transparency  [FLG: OFF]          : Invert the transparency logic (default: false)
   --alpha-constant-max   [INT: 255]          : Maximum value for constant alpha (default: 255)
   --alpha-quadratic-coef [FLT: 16.00]        : Coefficient for quadratic alpha calculation (default: 16.0)
   --alpha-linear-coef    [FLT: 1.00]         : Coefficient for linear alpha calculation (default: 1.0)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
