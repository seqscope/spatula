# spatula hist

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

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