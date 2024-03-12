# spatula build-sbcds

## Summary 

`spatula build-sbcds` creates a spatial barcode dictionary based from an Illumina FASTQ file. Here is a summary of the main features:

* Input: Takes a single-ended Illumina-sequenced FASTQ file containing the original read names as  input.
* Output: Produces a map of spatial barcodes (i.e. a mapping between barcode sequences and spatial coordinates), separated by each tile. It also generates a manifest file that summarizes the output for each tile. 


An example usage of the tool is as follows:

```sh
spatula build-sbcds --fq1 /path/to/input_file.R1.fastq.gz --format DraI32 \
                    --out /path/to/output/dir/
```

See below for a more detailed usage description.

## Required options
* `--fq` : The path to the (single-ended) input FASTQ file. The read name must include the spatial coordinate information. (i.e. SRA-processed read names must not be used)
* `--format` : Predefined spatial barcode format. It specifies the expected patterns of spatial barcodes, and whether reverse complement is required or not. A limited number of barcode patterns, including experimental ones, are accepted. The most popular choices are:
    -  `DraI32` : 32-bp barcode format from the Seq-Scope paper. The expected pattern in IUPAC codes is `NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN`; reverse complement is performed.
    - `HDMI20` : 20-bp barcode format from the original Seq-Scope paper. This is intended for a small array area; reverse complement is performed.
    - `DraI31` : 31-bp barcode format, almost identical to `DraI32`. This format is useful for one of the example sequence data we provided, which sequenced only the first 31-bp of `DraI32` format. The expected pattern in IUPAC codes is `NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBN`; reverse complement is performed.
* `--out` : Output directory that stores the spatial barcode dictionary per tile. See [Expected Output](#expected-output) for more details.

## Additional Options
* `--platform` : Platform on which the FASTQ file is generated. This specifies the expected format of the read name. The key fields extracted are `LANE`, `TILE`, `X`, `Y` coordinates. Currently, the only officially supported format is `Illumina` (Other experimental formats are not documented here and will not be supported).
    - `Illumina` : Illumina read name should have the following format: `[INSTRUMENT_ID]:[RUN_ID]:[FLOWCELL_ID]:[LANE]:[TILE]:[X]:[Y]`. 
* `--force-lane` : (For experimental platforms only) Force the lane number to a certain value; this is useful for certain platforms where the read name does not contain the lane information. 

## Expected Output

In the output directory `[outdir]`, the following files will be created.

* `[outdir]/manifest.tsv` contains the list of output tsv files. Each line contains the following information:
    - `id` : ID of the tile as `[lane]_[tile]`
    - `filepath` : Output file name associated with the tile `[lane]_[tile].sbcds.sorted.tsv.gz`. Note that the output is an unsorted tsv file, which must be sorted externally.
    - `barcodes` : The total number of spatial barcodes found in the tile.
    - `matches` : The number of barcodes matching the expected format.
    - `mismatches` : The number of barcodes that do not match the expected format.
    - `xmin` : The minimum value of the x coordinate per tile
    - `xmax` : The maximum value of the x coordinate per tile
    - `ymin` : The minimum value of the y coordinate per tile.
    - `ymax` : The maximum value of the y coordinate per tile.
* For each tile, `[outdir]/[lane]_[tile].sbcds.unsorted.tsv` contains the barcode sequences and their spatial coordinates in tsv format. Note that the output is an unsorted tsv file, which must be sorted alphabetically and compressed into `[lane]_[tile].sbcds.sorted.tsv.gz` after `build-sbcds` has finished. Each column in the tsv file contains the following information:
    1. Spatial barcode sequences (in reverse complement if the specified format requires it).
    2. Lane 
    3. Tile
    4. x-coordinate
    5. y-coordinate
    6. Number of bases that do not match the expected pattern defined by the format (0 is a perfect match).


## Full Usage 

The full usage of `spatula build-sbcds` can be viewed with the `--help` option:

```
$ ./spatula build-sbcds --help
[./spatula build-sbcds] -- Create spatial barcode dictionary based from HDMI FASTQ arrays

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --fq         [STR: ]             : Input FASTQ file that contains spatial barcode in the readname
   --format     [STR: ]             : Format of the spatial barcodes (e.g. DraI32, HDMI20, T7-30, etc)
   --platform   [STR: Illumina]     : Platform of the sequencing data to determine the rule to parse the readnames (e.g. Illumina)
   --force-lane [INT: 0]            : Force a lane number. Required a positive value for Salus/SalusGlobal platforms

== Output Options ==
   --out        [STR: ]             : Output directory name


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```