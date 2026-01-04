# spatula join-pixel-decode

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

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