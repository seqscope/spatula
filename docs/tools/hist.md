# spatula hist

## Summary 

`spatula hist` generates a text-based histogram from a specified numeric column of an input TSV file. It automatically or manually determines bin widths and calculates the frequency of data points within each bin. The tool can also report fractions, cumulative counts, and median values for each bin.

## Required Options

* `--file STR` : Input TSV file. Use `-` for stdin. (Default: `-`)
* `--out STR` : Output TSV file to show the results. Use `-` for stdout. (Default: `-`)

## Additional Options

* `--column INT` : 1-based index of the column to select for histogram generation. (Default: 1)
* `--bin-width FLT` : Fixed width of each bin.
* `--num-bins INT` : Number of bins to divide the data range into. (Default: 10). Note that `--bin-width` and `--num-bins` cannot be specified simultaneously.
* `--show-fraction` : If enabled, adds a `frac` column showing the fraction of total data points in each bin.
* `--show-cumulative` : If enabled, adds `cumul` (cumulative count) and optionally `fcumul` (cumulative fraction) columns.
* `--show-median` : If enabled, adds a `median` column showing the median value of data points within each bin.
* `--batch-size INT` : Size of the initial batch of records used to estimate the minimum and maximum values for determining bin width. (Default: 1,000,000)

## Expected Output

The output is a TSV file with the following columns:

* `from`: The lower bound of the bin interval.
* `to`: The upper bound of the bin interval.
* `median`: (Optional) The median value of the data points in the bin (only with `--show-median`).
* `count`: The number of data points falling into the bin.
* `frac`: (Optional) The fraction of total data points in the bin (only with `--show-fraction`).
* `cumul`: (Optional) The cumulative count of data points up to the current bin (only with `--show-cumulative`).
* `fcumul`: (Optional) The cumulative fraction of data points up to the current bin (only with `--show-cumulative` and `--show-fraction`).

## Full Usage 

The full usage of `spatula hist` can be viewed with the `--help` option:

```
$ ./spatula hist --help         
[./spatula hist] -- Create a text-based histogram

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --file            [STR: -]            : Input file. Use - for stdin
   --column          [INT: 1]            : 1-based index of the column to select

== Output options ==
   --out             [STR: -]            : Output tsv file to show the results. Use - for stdout
   --bin-width       [FLT: 0.00]         : Width of the each bin
   --num-bins        [INT: 0]            : Number of bins (10 by default)
   --show-fraction   [FLG: OFF]          : Show the fraction of the total
   --show-cumulative [FLG: OFF]          : Show the cumulative fraction
   --show-median     [FLG: OFF]          : Show the median value for each interval

== Other settings ==
   --batch-size      [INT: 1000000]      : Size of initial batch to determine the bin width


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```