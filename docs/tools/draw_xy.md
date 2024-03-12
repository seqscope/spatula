# spatula draw-xy

## Summary 

`spatula draw-xy` is a simple tool that draws an image of 2D points provided as an input. 

**IMPORTANT** [ImageMagick](https://imagemagick.org/script/download.php) must be installed to use this tool.

A typical use case is as follows:

* Input: TSV file (or can be streamed from STDIN), that contains X and Y coordinates of individual points to be plotted. 
* Output: Produces a grayscale 2D image that plots all the input points based on the specified scaling factor and intensity.

A typical example is as follows:

```sh
spatula draw-xy --tsv /path/to/input.tsv.gz \
                --out /path/to/output.png \
                --coord-per-pixel 1000 \
                --icol-x 3 \
                --icol-y 4 \
                --intensity-per-obs 50
```

See below for a more detailed usage description.

## Required options
* `--tsv` : The input TSV file that contains spatial coordinates of the points. `/dev/stdin` can be used to stream the input from STDIN.
* `--out` : The output filename of the image. Currently, `.png` is supported. 

## Additional Options

* `--width` : The width of the output image. If not specified, the width will be determined based on the minimum and maximum coordinates of the input points.
* `--height` : The height of the output image. If not specified, the height will be determined based on the minimum and maximum coordinates of the input points.
* `--coord-per-pixel` : The number of coordinates to be collapsed into a pixel as a factor to divide the input coordinate with. The default is 1.0.
* `--intensity-per-obs` : The intensity (max 255) of individual points to contribute to the pixel. The default is 1.0.
* `--icol-x` : The (0-based) column index of the X coordinate in the input TSV file. The default is 0.
* `--icol-y` : The (0-based) column index of the Y coordinate in the input TSV file. The default is 1.

## Expected Output

The output `[out]` will be created as a PNG file containing the image of the input points.

## Full Usage 

The full usage of `spatula draw-xy` can be viewed with the `--help` option:

```
$ ./spatula draw-xy --help           
[./spatula draw-xy] -- Draw the image of points in 2D space

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --tsv               [STR: ]             : tsv file to draw the x-y coordinates. /dev/stdin for stdin
   --icol-x            [INT: 0]            : 0-based index of the column for x
   --icol-y            [INT: 1]            : 0-based index of the column for y

== Settings ==
   --width             [INT: 0]            : Width of the image
   --height            [INT: 0]            : Height of the image
   --coord-per-pixel   [FLT: 1.00]         : Number of coordinate units per pixel
   --intensity-per-obs [INT: 1]            : Intensity per pixel per observation

== Output Options ==
   --out               [STR: ]             : Output file name


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```