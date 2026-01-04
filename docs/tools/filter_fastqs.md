# spatula filter-fastqs

## Summary 

TBA

## Required options

TBA

## Additional Options

TBA 

## Expected Output

TBA

## Full Usage 

The full usage of `spatula filter-fastqs` can be viewed with the `--help` option:

```
$ ./spatula filter-fastqs --help         
[./spatula filter-fastqs] -- Filter FASTQs based on a pattern

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --fq1          [STR: ]             : FASTQ file for read 1
   --fq2          [STR: ]             : FASTQ file for read 2

== Filtering options ==
   --pat1         [V_STR: ]           : IUPAC pattern to match for Read 1 (1+ match required)
   --pat2         [V_STR: ]           : IUPAC pattern to match for Read 2 (1+ match required)
   --min-mismatch [INT: 0]            : Minimum number of mismatches allowed per pattern
   --remove       [FLG: OFF]          : Remove the matching sequence from the FASTQ file instead of keeping it

== Output Options ==
   --out1         [STR: ]             : Output FASTQ file for read1
   --out2         [STR: ]             : Output FASTQ file for read2


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```