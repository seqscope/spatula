# spatula draw-sge

## Summary 

`spatula draw-sge` is a tool to visualize the spatial gene expression based on specific combination of colors and gene set. 

**IMPORTANT** [ImageMagick](https://imagemagick.org/script/download.php) must be installed to use this tool.

A typical use case is as follows:

* Input: Takes a SGE matrix from `dge2sge` and the list of genes (with designated colors) to visualize the spatial gene expression.
* Output: Produces an 2D image that plots the spatial distribution of expression of selected genes.

A typical example is as follows:

```
spatula draw-sge  --manifest /path/to/combine/sbcds/output/dir/manifest.tsv \
                  --sge      /path/to/dge2sge/output/dir \
                  --color-gene 320000:_all_:1 \ 
                  --color-gene 003200:Glul,Cyp2e1:1 \
                  --color-list 000064:/path/to/custom/gene_list.tsv \
                  --out      /path/to/output/image.png
                  --auto-adjust                
```

See below for a more detailed usage description.

## Required options

* `--sge` : The directory containing SGE matrix, typically from `dge2sge` command.
* `--out` : The output filename of the image. Currently, `.png` is supported. 

## Options to Specify Genes to Visualize

* `--color-gene` : (Allows multiple uses) A string formatted as [RGB_Hex_Code]:[gene1],[gene2],...:[idx]. 
    - The `[RGB_Hex_Code]` is the RGB hex color code (`RRGGBB` format) for each observation of the specified genes.
    - The `[gene1],[gene2],...` is the list of genes to color with the specified color. Either gene ID (e.g. `ENSMUSG00000026473`) or gene symbol (e.g. `Glul`) may be used. `_all_` can be used to color all genes with the specified color.
    - The `[idx]` is optional field to specify the index of the gene expression in the SGE matrix. (e.g. 1: Gene, 2: GeneFull, 3: Spliced, 4: Unspliced, 5: Velocyto). The default is 1.
    - For example, `--color-gene 003200:Glul,Cyp2e1:1` will color *Glul* and *Cyp2e1* genes with green (increasing intensity by 50 for each observation).
* `--color-regex` : (Allows multiple uses) A string formatted as [RGB_Hex_Code]:[regex_pattern]:[idx]. 
    - The `[RGB_Hex_Code]` is the RGB hex color code (`RRGGBB` format) for each observation of the specified genes.
    - The `[regex_pattern` is the regulat expression pattern to specify the list of genes. Here are some examples:
      - You can specify multiple genes like `^(Glul|Cyp2e1)$` 
      - You can specify genes with a common prefix like `^mt-.*$` (for mitochondrial genes)
    - The `[idx]` is optional field to specify the index of the gene expression in the SGE matrix. (e.g. 1: Gene, 2: GeneFull, 3: Spliced, 4: Unspliced, 5: Velocyto). The default is 1.
    - For example, `--color-gene 003200:^(Glul|Cyp2e1):1` will color *Glul* and *Cyp2e1* genes with green (increasing intensity by 50 for each observation).
* `--color-list` : (Allows multiple uses) A string formatted as [RGB_Hex_Code]:[path_to_gene_list]:[default_idx].
    - The `[RGB_Hex_Code]` is the RGB hex color code (`RRGGBB` format) for each observation of the specified genes.
    - The `[path_to_gene_list]` is the path to the gene list file. The file should contain a list of gene names of gene IDs (without headers)
    - The `[default_idx]` is optional field to specify the index of the gene expression in the SGE matrix. (e.g. 1: Gene, 2: GeneFull, 3: Spliced, 4: Unspliced, 5: Velocyto). This is only in effect when `idx` column is not present in the `[path_to_gene_list]` file.
    - If `[default_idx]` is absent, it will be defaulted to 1 (Gene)
    - For example, `--color-list 000064:/path/to/custom/gene_list.tsv` will color genes in the `/path/to/custom/gene_list.tsv` file with blue, with `default_idx=1`, increasing intensity by 100 for each observation.
    - Here is an example content of the gene list file `mt-genes.tsv`, which contains gene ID for mitochondrial genes in mouse:
```
mt-Nd1
mt-Nd2
mt-Co1
mt-Co2
mt-Atp8
mt-Atp6
mt-Co3
mt-Nd3
mt-Nd4l
mt-Nd4
mt-Nd5
mt-Nd6
mt-Cytb
```

## Additional Options

* `--manifest` : The `manifest.tsv` file from the `combine-sbcds` file that contains the summary of the spatial coordinate of a Seq-Scope Chip. It must contain `xmin`, `xmax`, `ymin`, and `ymax` columns. This option is REQUIRED unless `--auto-adjust` is used.
* `--auto-adjust` : Automatically adjust the intensity of the color based on the maximum count. The default is OFF.
* `--adjust-quantile` : The quantile of pixel to use for auto-adjustment among non-zero pixels. The default is 0.99.
* `--coord-per-pixel` : The number of coordinates to be collapsed into a pixel as a factor to divide the input coordinate with. The default is 1000.0.
* `--bcd` : The barcode file name in the SGE directory. The default is `barcodes.tsv.gz`.
* `--ftr` : The feature file name in the SGE directory. The default is `features.tsv.gz`.
* `--mtx` : The matrix file name in the SGE directory. The default is `matrix.mtx.gz`.

## Expected Output

The output `[out]` will be created as a PNG file containing the image of the input points.

## Full Usage 

The full usage of `spatula draw-sge` can be viewed with the `--help` option:

```
$ ./spatula draw-sge --help
[./spatula draw-sge] -- Draw the image of spatial gene expression (SGE) data

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input files ==
   --minmax          [STR: ]             : Bounding box information. Expects xmin/xmax/ymin/ymax (tall or wide format)
   --sge             [STR: ]             : SGE directory
   --bcd             [STR: barcodes.tsv.gz] : Barcode file name
   --ftr             [STR: features.tsv.gz] : Feature file name
   --mtx             [STR: matrix.mtx.gz] : Matrix file name

== Genes to visualize ==
   --color-gene      [V_STR: ]           : [color_code]:[gene1],[gene2],... as a visualization unit. Adding :[idx] at the end is optional
   --color-regex     [V_STR: ]           : [color_code]:[regex](:[idx]) as a visualization unit. [regex] is a regulat expression. Adding :[idx] at the end is optional
   --color-list      [V_STR: ]           : [color_code]:[list_file](:[idx]) as a visualization unit

== Output options ==
   --coord-per-pixel [FLT: 1000.00]      : Number of coordinate units per pixel
   --auto-adjust     [FLG: OFF]          : Automatically adjust the intensity of the color based on the maximum count
   --adjust-quantile [FLT: 0.99]         : Quantile of pixel to use for auto-adjustment among non-zero pixels

== Output Options ==
   --out             [STR: ]             : Output file name

NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```