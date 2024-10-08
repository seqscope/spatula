# spatula convert-sge

## Summary 

`spatula convert-sge` convert digital gene expression (DGE) data following the [10x Genomics Cell Ranger Feature Barcode Matrix format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices) into a generic spatial gene expression TSV file format that can be analyzed by [FICTURE](https://github.com/seqscope/ficture). It can also convert spatial gene expression (SGE) data generated by the [NovaScope](https://github.com/seqscope/novascope) pipeline into the same generic spatial gene expression TSV file format.

Here is a summary of the main features:

* Input: Takes (a) `barcode.tsv.gz` file containing spatial barcodes, (b) `features.tsv.gz` file containing gene names, and (c) `matrix.mtx.gz` file containing the sparsely encoded gene expression matrix following the [MatrixMarket](https://math.nist.gov/MatrixMarket/formats.html) format. If the spatial barcodes do not contain the spatial coordinates (e.g. 10x Visium HD), additional barcode position files must be provided separately in CSV format. 
* Output: A TSV file containing the spatial gene expression data in the format of `gene`, `barcode`, `x`, `y`, and `count`. The spatial coordinates are converted to micrometer scale. Note that the output file is NOT sorted by the spatial coordinates.

## Usage

### Converting Seq-Scope/NovaScope Feature Barcode Matrix

The standard output from the [SeqScope Protocol](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11014489/) using the [NovaScope](https://github.com/seqscope/novascope) pipeline generates the following set of output files, assuming `${SGEDIR}` is the directory containing the SGE output from NovaScope. 

* `${SGEDIR}/barcodes.tsv.gz` : The barcode file containing the spatial barcodes.
* `${SGEDIR}/features.tsv.gz` : The feature file containing the gene names.
* `${SGEDIR}/matrix.mtx.gz` : The matrix file containing the gene expression matrix in the MatrixMarket format.

Note that the `barcodes.tsv.gz` file contains spatial coordinates in their output. For example, in NovaScope 1.0, the X/Y coordinates of each spatial barcode is located in the 6th and 7th columns in nanometer scale, respectively.

You may obtain the example output of the NovaScope output at our [Zenodo repository](https://zenodo.org/records/12773392).

```bash
## Show the first few lines of the barcode file
$ gzip -cd ${SGEDIR}/barcodes.tsv.gz | head

AAAAAAAATAGTTCTGCTAGCTGGTAAGCT  1       1       1       1       7124822 2906042 1,1,1,0,0
AAAAAAAGTGATCAGAGGTGATATTATGCT  2       2       1       1       7382402 2717083 3,3,3,0,0
AAAAAAAGTTCGCACTATACGAACAGGGAT  3       3       1       1       8634969 2839091 1,1,1,0,0
AAAAAAGGTACCCGCAGTGCGGACAAACGA  4       4       1       1       4994894 2666443 1,1,1,0,0
AAAAACGCTCCCCGTTGTGAATGGGGTCGT  5       5       1       1       8575635 3991700 1,1,1,0,0
AAAAACTCTCAGAAGAGAAAGTAATAGTCG  6       6       1       1       6747424 2860446 1,0,0,0,0
AAAAAGAGAACCACAGGTAATCCACCTACA  7       7       1       1       5771532 3028677 5,4,3,0,0
AAAAAGAGGTGAGGGTCGCCTGCATATTAG  8       8       1       1       6614008 3251930 2,2,2,0,0
AAAAAGGGGTCTAGAGGAGACAATGAAGTG  9       10      1       1       4634991 4060530 2,2,2,0,0
AAAAATTATGGACGACCTACTTCTCGGTGG  10      15      1       1       5358763 4389131 1,1,1,0,0
```

To convert Seq-Scope/NovaScope output to a FICTURE-compatible format, the following two steps are needed.

#### 1. Run `spatula convert-sge` to transform the input files into a generic TSV file.

Note that you need to install the `spatula` tool before running the following command.
Because the spatial coordinates are in nanometer scale, the `--units-per-um` parameter should be set to 1000 to convert the coordinates to micrometer scale.

```bash
## Run spatula convert-sge, assuming that output_dir is the output directory
$ spatula convert-sge \
    --in-sge ${SGEDIR} \
    --units-per-um 1000 \
    --colnames-count Count \
    --out output_dir 
## Uncomment the following line if you want to exclude commonly ignored features
#   --exclude-feature-regex '^(BLANK|NegCon|NegPrb|mt-|MT-|Gm\d+$$)'
## The following parameters follow the default values (unnecessary to specify)
#   --sge-bcd barcode.tsv.gz \
#   --sge-ftr features.tsv.gz \
#   --sge-mtx matrix.mtx.gz \
#   --icols-mtx 1 \
#   --icol-bcd-barcode 1 \
#   --icol-bcd-x 6 \
#   --icol-bcd-y 7 \
#   --icol-ftr-id 1 \
#   --icol-ftr-name 2 \
#   --colname-x X \
#   --colname-y Y \
#   --colname-feature-name gene 
```

Make sure that the output file is generated in the expected format. 

```
$ gzip -cd out_dir/transcripts.unsorted.tsv.gz | head

X       Y       gene    Count
7124.82 2906.04 Pmpcb   1
7382.40 2717.08 Acin1   1
7382.40 2717.08 Slc35b1 1
8634.97 2839.09 Marchf7 1
4994.89 2666.44 Exosc9  1
8575.64 3991.70 Gpld1   1
6747.42 2860.45 Serpina1d       1
5771.53 3028.68 Ccnl2   1
5771.53 3028.68 Ugt2b34 1
```

#### 2. Sort the generic TSV file by a spatial coordinate.

You may sort the output file by the spatial coordinates using the `sort` command. For example, to sort the file by the first column (X-coordinate), you may run the following command.

```bash
## Sort the unsorted output file by the X-coordinate
## Note that this step will take a while for a large file
(gzip -cd output_dir/transcripts.unsorted.tsv.gz \
    | head -1; gzip -cd output_dir/transcripts.unsorted.tsv.gz \
    | tail -n +2 | sort -S 1G -gk1) \
    | gzip -c > output_dir/transcripts.sorted.tsv.gz
```


### Converting 10x Genomics Visium HD Feature Barcode Matrix

The standard output from the 10x Genomics Visium HD pipeline contains many files. Typically, the following files are used for the conversion, assuming `${DATADIR}` is the directory containing the Visium HD data. 

* `${DATADIR}/square_002um/raw_feature_bc_matrix/barcodes.tsv.gz` : The barcode file containing the spatial barcodes.
* `${DATADIR}/square_002um/raw_feature_bc_matrix/features.tsv.gz` : The feature file containing the gene names and IDs.
* `${DATADIR}/square_002um/raw_feature_bc_matrix/matrix.mtx.gz` : The matrix file containing the gene expression matrix in the MatrixMarket format.
* `${DATADIR}/square_002um/spatial/tissue_positions.parquet` : The barcode position file containing the spatial coordinates of the barcodes in parquet format.
* `${DATADIR}/square_002um/spatial/scalefactors_json.json` : The scale factor file containing the pixel-to-micrometer conversion factor `microns_per_pixel`.

In this example, we used the FFPE mouse brain data
downloaded at [10x Genomics website](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he).

Instead of using `raw_feature_bc_matrix` directory `filtered_feature_bc_matrix` directory may be used if the filtered gene expression matrix is preferred. 

To convert Visium HD data to a FICTURE-compatible format, the following four steps are needed.

#### 1. Obtain `units-per-um` parameter from the `scalefactors_json.json` file.

For example, if you may find the `microns_per_pixel` value as follows:

```bash
## Find microns_per_pixel value from the scalefactors_json.json file
$ grep -w microns_per_pixel ${DATADIR}/square_002um/spatial/scalefactors_json.json 
    "microns_per_pixel": 0.2738242950835738,
```

In this case, the `units-per-um` value is `1/0.2738242950835738 = 3.652`.

#### 2. Convert `tissue_positions.parquet` into a CSV file.

The `tissue_positions.parquet` file can be converted into a CSV file using the `parquet-tools` command. 
Please install [parquet-tools](https://pypi.org/project/parquet-tools/) before running the following command.

```sh
## Convert tissue_positions.parquet into a CSV file
$ parquet-tools csv ${DATADIR}/square_002um/spatial/tissue_positions.parquet \
    | gzip -c > ${DATADIR}/square_002um/spatial/tissue_positions.csv.gz

## Check the first few lines of the CSV file
$ gzip -cd 10x_visium_hd/raw/tissue_positions.csv.gz | head

barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
s_002um_00000_00000-1,0,0,0,44.1895274751605,21030.743287971593
s_002um_00000_00001-1,0,0,1,44.30210038790203,21023.440348914042
s_002um_00000_00002-1,0,0,2,44.4146732713411,21016.137411757423
s_002um_00000_00003-1,0,0,3,44.52724612547773,21008.83447650175
s_002um_00000_00004-1,0,0,4,44.639818950311934,21001.531543147008
s_002um_00000_00005-1,0,0,5,44.7523917458437,20994.228611693197
s_002um_00000_00006-1,0,0,6,44.86496451207307,20986.92568214033
s_002um_00000_00007-1,0,0,7,44.977537249000044,20979.622754488388
s_002um_00000_00008-1,0,0,8,45.09010995662463,20972.319828737378
```

#### 3. Run `spatula convert-sge` to transform the input files into a generic TSV file.

Note that you need to install the `spatula` tool before running the following command.

```bash
## Run spatula convert-sge, assuming that output_dir is the output directory
$ spatula convert-sge \
    --in-sge ${DATADIR}/square_002um/raw_feature_bc_matrix \
    --pos ${DATADIR}/square_002um/spatial/tissue_positions.csv.gz \
    --units-per-um 3.652 \
    --colnames-count Count \
    --out output_dir 
## Uncomment the following line if you want to exclude commonly ignored features
#   --exclude-feature-regex '^(BLANK|NegCon|NegPrb|mt-|MT-|Gm\d+$$)'
## The following parameters follow the default values (unnecessary to specify)
#   --sge-bcd barcode.tsv.gz \
#   --sge-ftr features.tsv.gz \
#   --sge-mtx matrix.mtx.gz \
#   --icols-mtx 1 \
#   --icol-bcd-barcode 1 \
#   --icol-ftr-id 1 \
#   --icol-ftr-name 2 \
#   --pos-colname-barcode barcode \
#   --pos-colname-x pxl_row_in_fullres \
#   --pos-colname-y pxl_col_in_fullres \
#   --pos-delim , \
#   --colname-x X \
#   --colname-y Y \
#   --colname-feature-name gene 
```

Make sure that the output file is generated in the expected format. 

```
$ gzip -cd output_dir/transcripts.unsorted.tsv.gz | head

X       Y       gene    Count
13.58   5662.74 Ctps2   1
14.04   5632.75 Rab28   1
14.19   5752.76 Omg     1
14.19   5752.76 Prnp    1
14.99   5700.77 Sorbs1  1
15.70   5654.77 Myoc    1
15.89   5642.78 Fgfr1op2        1
16.53   5600.78 Hmgcr   1
18.32   5484.80 Stmn1   1
```

#### 4. Sort the generic TSV file by a spatial coordinate.

You may sort the output file by the spatial coordinates using the `sort` command. For example, to sort the file by the first column (X-coordinate), you may run the following command.

```bash
## Sort the unsorted output file by the X-coordinate
## Note that this step will take a while for a large file
(gzip -cd output_dir/transcripts.unsorted.tsv.gz \
    | head -1; gzip -cd output_dir/transcripts.unsorted.tsv.gz \
    | tail -n +2 | sort -S 1G -gk1) \
    | gzip -c > output_dir/transcripts.sorted.tsv.gz
```

## Command line arguments
### Key options

* `--in-sge` : The path to the directory containing spatial barcode files.
* `--out-tsv` : The path to the directory to store the output files.
* `--icols-mtx` : Comma-separated 1-based column indices use as the count in the `matrix.mtx.gz` file.
* `--colnames-count` : Comma-separate column names for Count for each column in the output file. The order should match to `--icols-mtx`.
* `--icol-bcd-barcode` : The 1-based column index of the barcode in the `barcodes.tsv.gz` file.
* `--icol-bcd-x` : 1-based column index of X coordinate in the `barcodes.tsv.gz` file (when `--pos` is not provided).
* `--icol-bcd-y` : 1-based column index of Y coordinate in the `barcodes.tsv.gz` file (when `--pos` is not provided).
* `--icol-ftr-id` : 1-based column index of the feature (gene) ID in the `features.tsv.gz` file.
* `--icol-ftr-name` : 1-based column index of the feature name (gene symbol) in the `features.tsv.gz` file.
* `--pos` : The path to the barcode position file (typically used for Visium HD input in CSV format).
* `--units-per-um` : The number of units (i.e. pixels) per micrometer for the spatial coordinates.

### Gene filtering options
* `--include-feature-list` ; A file containing a list of input genes to be included (feature name of IDs)
* `--exclude-feature-list` ; A file containing a list of input genes to be excluded (feature name of IDs)
* `--include-feature-substr` : A substring of feature/gene names to be included
* `--exclude-feature-substr` : A substring of feature/gene names to be excluded
* `--include-feature-regex` : A regular expression of feature/gene names to be included
* `--exclude-feature-regex` : A regular expression of feature/gene names to be excluded

### Additional options
* `--out-sge` : The path to the directory to store the output files in the barcode-feature-matrix format (*still experimental*).
* `--sge-bcd` : The name of the barcode file in the barcode-feature-matrix format. Default is `barcodes.tsv.gz`.
* `--sge-ftr` : The name of the feature file in the barcode-feature-matrix format. Default is `features.tsv.gz`.
* `--sge-mtx` : The name of the matrix file in the barcode-feature-matrix format. Default is `matrix.mtx.gz`.
* `--tsv-mtx` : The name of the matrix file in the generic TSV format. Default is `transcripts.unsorted.tsv.gz`.
* `--tsv-ftr` : The name of the feature file in the generic TSV format. Default is `features.tsv.gz`.
* `--tsv-minmax` : The name of the minmax coordinate file in the generic TSV format. Default is `minmax.tsv`.
* `--pos-colname-barcode` : The column name of the barcode in the barcode position file.
* `--pos-colname-x` : The column name of the X coordinate in the barcode position file.
* `--pos-colname-y` : The column name of the Y coordinate in the barcode position file.
* `--pos-delim` : The delimiter used in the barcode position file. Default is `,`.
* `--precision-um` : Output precision below the decimal point.
* `--print-feature-id` : Turn on to print feature ID in addition to feature name.
* `--allow-duplicate-gene-names` : Turn on to allow duplicate gene names in the output file.
* `--colname-feature-name` : Column name for feature/gene name in the ouput file (default is `gene`).
* `--colname-feature-id` : Column name for feature/gene ID (default is `gene_id`).
* `--colname-x` : Column name for X-coordinate in the output file (default is `X`).
* `--colname-y` : Column name for Y-coordinate in the output file (default is `Y`).

## Expected Output

In the output directory `[outdir]`, the following files will be created.

- `[outdir]/transcript.unsorted.tsv.gz` : A generic TSV file that contains individual spatial coordinates and transcript counts. Each column in the tsv file contains the following information. The column names may vary depending on the input arguments.
    - `X` : The X-coordinate of the spatial barcode / molecule in micrometer.
    - `Y` : The Y-coordinate of the spatial barcode / molecule in micrometer.
    - `gene` : The gene name of the transcript.
    - `Count` : The number of transcripts detected in the spatial barcode. May be named as `gn`, `gt` or other names depending on the `--colnames-count` argument. 
- `[outdir]/feature.tsv.gz` : A TSV containing the list of gene names, possibly gene IDs, and overall counts per each gene.
- `[outdir]/minmax.tsv` : A TSV containing the minimum and maximum coordinates of the spatial barcodes.
    - `xmin` : The minimum X-coordinate of the spatial barcodes in micrometer.
    - `xmax` : The maximum X-coordinate of the spatial barcodes in micrometer.
    - `ymin` : The minimum Y-coordinate of the spatial barcodes in micrometer.
    - `ymax` : The maximum Y-coordinate of the spatial barcodes in micrometer.

## Full Usage 

The full usage of the software tool is as follows:

```sh
$ ./spatula convert-sge --help
[./spatula convert-sge] -- Convert SGE into plain various formats

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Key Input/Output Options ==
   --in-sge                     [STR: ]             : Input SGE directory
   --out-sge                    [STR: ]             : Path of SGE output directory
   --out-tsv                    [STR: ]             : Path of TSV output directory

== File names for SGE/TSV input/output ==
   --sge-bcd                    [STR: barcodes.tsv.gz] : Barcode file name in SGE directory
   --sge-ftr                    [STR: features.tsv.gz] : Feature file name in SGE directory
   --sge-mtx                    [STR: matrix.mtx.gz] : Matrix file name in SGE directory
   --tsv-mtx                    [STR: transcripts.unsorted.tsv.gz] : Transcript file name in TSV directory
   --tsv-ftr                    [STR: features.tsv.gz] : Feature file name in TSV directory
   --tsv-minmax                 [STR: minmax.tsv]   : Minmax file name in TSV directory

== Expected column index in SGE input ==
   --icols-mtx                  [STR: 1,2,3,4,5]    : Comma-separated 1-based column indices use as the count
   --icol-bcd-barcode           [INT: 1]            : 1-based column index of barcode in the barcode file
   --icol-bcd-x                 [INT: 6]            : 1-based column index of x coordinate in the barcode file
   --icol-bcd-y                 [INT: 7]            : 1-based column index of y coordinate in the barcode file
   --icol-ftr-id                [INT: 1]            : 1-based column index of feature ID in the barcode file
   --icol-ftr-name              [INT: 2]            : 1-based column index of feature name in the barcode file

== Additional Barcode Position File ==
   --pos                        [STR: ]             : Position file name that contains separate X and Y coordinates
   --pos-colname-barcode        [STR: barcode]      : Column name for barcode in the position file
   --pos-colname-x              [STR: pxl_row_in_fullres] : Column name for X-axis in the position file
   --pos-colname-y              [STR: pxl_col_in_fullres] : Column name for Y-axis in the position file
   --pos-delim                  [STR: ,]            : Delimiter for the position file (default: ",")

== Input Filtering Options ==
   --include-feature-list       [STR: ]             : A file containing a list of input genes to be included (feature name of IDs)
   --exclude-feature-list       [STR: ]             : A file containing a list of input genes to be excluded (feature name of IDs)
   --include-feature-substr     [STR: ]             : A substring of feature/gene names to be included
   --exclude-feature-substr     [STR: ]             : A substring of feature/gene names to be excluded
   --include-feature-regex      [STR: ]             : A regex pattern of feature/gene names to be included
   --exclude-feature-regex      [STR: ]             : A regex pattern of feature/gene names to be excluded

== Key Output Options ==
   --units-per-um               [FLT: 1.00]         : Coordinate unit per um (conversion factor)
   --precision-um               [INT: 3]            : Output precision below the decimal point

== Auxilary Output Options for TSV output ==
   --print-feature-id           [FLG: OFF]          : Print feature ID in output file
   --allow-duplicate-gene-names [FLG: OFF]          : Allow duplicate gene names in the output file
   --colname-feature-name       [STR: gene]         : Column name for feature/gene name
   --colname-feature-id         [STR: gene_id]      : Column name for feature/gene ID
   --colnames-count             [STR: gn,gt,spl,unspl,ambig] : Comma-separate column names for Count
   --colname-x                  [STR: X]            : Column name for X
   --colname-y                  [STR: Y]            : Column name for Y


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```