# spatula custom-demux-fastq

## Summary 

`spatula custom-demux-fastq` takes FASTQ files of sequenced reads and indexes to demultiplex the reads by sample. Compared to the default FASTQ demultiplexing tool provided by Illumina, `custom-demux-fastq` allows users to specify the matching criteria in a more flexible way. Specifically, the software performs demultiplexing based on the following criteria:

* Unlike the default Illumina demultiplexing tool, `custom-demux-fastq` does not consider `N` as mismatches (by default).
* The software allows users to specify the number of mismatches allowed to be considered as the match.
* It also requires that the 2nd-best matching index sequences have larger hamming distance with the index sequences with a specific margin, compared to the best matching index sequences.

Note that this tool is NOT limited to demultiplexing spatial transcriptomics sequence data, and can be used for a general tool for demultiplexing FASTQ files as well.

Also note that this tool uses multi-process parallelization to speed up the compression. It is recommended to run this tool with at least 5-10 CPUs. `gzip` or `pigz` must be installed in the system to compress the output files.

Here is a summary of a typical use case:

* Input: Takes FASTQ files that contain (a) Read 1, (b) Read 2 (optional), (c) Index 1, (d) Index 2 (optional) sequences.
* Output: Produces multiple sets of FASTQs separated by the sample identity.

A typical example is as follows:

```sh
spatula custom-demux-fastq --R1 /path/to/input.R1.fastq.gz \
                           --R2 /path/to/input.R2.fastq.gz \
                           --I1 /path/to/input.I1.fastq.gz \
                           --sample /path/to/sample.index.tsv \
                           --out /path/to/output_prefix
```

See below for a more detailed usage description.

## Detailed Usage Description

There are two expected use cases of `spatula custom_demux_fastq`:

### Demultiplexing from BCL files

If you have access to the full BCL file, instead of demultiplexing individual samples, you may run `bcl2fastq` tool without performing demultiplexing, but creating index sequences as follows: 

```bash
## ${bcldir} : directory containing BCL files
## ${outdir} : directory of output FASTQ files
bcl2fastq -R ${bcldir} -o ${outdir} --create-fastq-for-index-reads
```

Then all reads in the FASTQ file will be written into "Undetermined" FASTQs. You can use `spatula custom_demux_fastq` as input to demultiplex FASTQ files instead. 

### Further demultiplexing already demultiplexed FASTQ files

If you already have demultiplexed FASTQ files, you will already have FASTQ files that are demultiplexed individual samples. We typically do NOT recommend running `spatula custom_demux_fastq` on the FASTQ files that are successfully demultiplexed, as the results will look very similar to the default `bcl2fastq` pipeline. 

However, if you have a substantial amount of "Undetermined" reads remaining, you may want use `spatula custom_demux_fastq` to further demultiplex the reads. Because Illumina's `bcl2fastq` pipeline typically performs demulitplexing in a conservative way, you may be able to rescue some of the reads with this tool.

Note that, if you have a very large number of samples (e.g. >10) demultiplexed in a single run, modifying the default parameter (e.g. using `--min-diff 1` and/or `--max-mismatch 1`) may be necessary to achieve more sensible results.


## Key options
* `--R1` : (Required) The path to the FASTQ file of Read 1.
* `--R2` : (Optional) The path to the FASTQ file of Read 2.
* `--I1` : (Required) The path to the FASTQ file of Index 1.
* `--I2` : (Optional) The path to the FASTQ file of Index 2.
* `--sample` : (Required) Sample index file in `tsv` format, with the following three columns:
     - `ID` : Sample ID. Must be unique.
     - `I1` : Index 1 sequence.
     - `I2` : Index 2 sequence (if available). 
* `--out` : (Required) The prefix of the output file ([out_prefix].[sample_ID][suffix]) paths. See [Expected Output](#expected-output) for more details.

## Additional Options

* `--cmd` : The command to run for compressing output files. Each file is compressed with the command specified. The default is `gzip -c`. Another popular choice could be `pigz -p 4 -c`, which is faster than `gzip -c` if `pigz` is installed.
* `--max-mismatch` : The maximum number of mismatches allowed to be considered as a match. The default is 2.
* `--min-diff` : Minimum difference in the hamming distance with the index sequence between the best and second best match. The default is 2.
* `--suffix-R1` : The suffix of Read 1 output files. The default is `.R1.fastq.gz`.
* `--suffix-R2` : The suffix of Read 2 output files. The default is `.R2.fastq.gz`.
* `--suffix-I1` : The suffix of Index 1 output files. The default is `.I1.fastq.gz`.
* `--suffix-I2` : The suffix of Index 2 output files. The default is `.I2.fastq.gz`.
* `--ambiguous` : The sample ID of the ambiguously demultiplexed reads. The default is `ambiguous`.

## Expected Output

With `[out_prefix]` as the prefix, and `[sample_ID]` be the ID of best-matching samples, the following files will be created:

* `[out_prefix].[sample_ID].R1.fastq.gz` : The FASTQ file of Read 1 for the sample.
* `[out_prefix].[sample_ID].R2.fastq.gz` : The FASTQ file of Read 2 for the sample (if available).
* `[out_prefix].[sample_ID].I1.fastq.gz` : The FASTQ file of Index 1 for the sample.
* `[out_prefix].[sample_ID].I2.fastq.gz` : The FASTQ file of Index 2 for the sample (if available).
* If the suffix is specified with `--suffix-R1`, `--suffix-R2`, `--suffix-I1`, `--suffix-I2`, the suffix will be replaced with the specified suffix.

## Full Usage 

The full usage of `spatula subset-sge` can be viewed with the `--help` option:

```
$ ./spatula custom-demux-fastq --help
[./spatula custom-demux-fastq] -- Demultiplex FASTQ files based in a customized manner

 Copyright (c) 2022-2024 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --R1                     [STR: ]             : FASTQ file for read 1
   --R2                     [STR: ]             : FASTQ file for read 2
   --I1                     [STR: ]             : FASTQ file for index 1
   --I2                     [STR: ]             : FASTQ file for index 2
   --sample                 [STR: ]             : Sample sheet file containing [ID] [I1] [I2]

== Settings ==
   --cmd                    [STR: gzip -c]      : Command to compress the output files
   --consider-N-as-mismatch [FLG: OFF]          : Consider N as mismatch (default: false)
   --max-mismatch           [INT: 2]            : Maximum number of mismatch allowed
   --min-diff               [INT: 2]            : Minimum difference between the best and second best match

== Output Options ==
   --out                    [STR: ]             : Output prefix
   --suffix-R1              [STR: .R1.fastq.gz] : Output suffix for read 1
   --suffix-R2              [STR: .R2.fastq.gz] : Output suffix for read 2
   --suffix-I1              [STR: .I1.fastq.gz] : Output suffix for index 1
   --suffix-I2              [STR: .I2.fastq.gz] : Output suffix for index 2
   --ambiguous              [STR: ambiguous]    : The keyword for ambiguous samples (default: ambiguous)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```