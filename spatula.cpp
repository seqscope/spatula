#include "spatula.h"
#include "qgenlib/commands.h"
#include "qgenlib/qgen_utils.h"

int32_t cmdMatchTag(int32_t argc, char **argv);
int32_t cmdSearchTag(int32_t argc, char **argv);
int32_t cmdMakeSpatialDGE(int32_t argc, char **argv);
int32_t cmdMatchSpatialBarcodes(int32_t argc, char **argv);
int32_t cmdMergeMatchedTags(int32_t argc, char **argv);
int32_t cmdBuildSpatialBarcodeDict(int32_t argc, char **argv);
int32_t cmdReformatFASTQs(int32_t argc, char **argv);
int32_t cmdFilterFASTQs(int32_t argc, char **argv);
int32_t cmdReformatFASTQsMT(int32_t argc, char **argv);
int32_t cmdDGE2SDGE(int32_t argc, char **argv);
int32_t cmdSGE2TSV(int32_t argc, char **argv);
int32_t cmdConvertSGE(int32_t argc, char **argv);
int32_t cmdSegmentSGE(int32_t argc, char **argv);
int32_t cmdDrawSGE(int32_t argc, char **argv);
int32_t cmdSubsetSGE(int32_t argc, char **argv);
int32_t cmdCombineSGE(int32_t argc, char **argv);
int32_t cmdEvalDupsSBCD(int32_t argc, char **argv);
int32_t cmdDrawXY(int32_t argc, char **argv);
int32_t cmdDraw3way(int32_t argc, char **argv);
int32_t cmdCombineSBCD(int32_t argc, char **argv);
int32_t cmdHist(int32_t argc, char **argv);
int32_t cmdFilterCommonBarcodes(int32_t argc, char **argv);
int32_t cmdStratifyFASTQByBarcodes(int32_t argc, char **argv);
int32_t cmdCustomDemuxFASTQ(int32_t argc, char **argv);
int32_t cmdJoinPixelTSV(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
{
  commandHelp.copyright_str = "Copyright (c) 2022-2024 by Hyun Min Kang";

  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
  LONG_COMMAND_GROUP("Tools for Analyzing Spatial Barcodes", NULL)
  LONG_COMMAND("build-sbcds", &cmdBuildSpatialBarcodeDict, "Create spatial barcode dictionary based from HDMI FASTQ arrays")
  LONG_COMMAND("combine-sbcds", &cmdCombineSBCD, "Combine multiple SBCD files")
  LONG_COMMAND("match-sbcds", &cmdMatchSpatialBarcodes, "Match the FASTQ file containing spatial barcodes with the spatial barcode dictionary")
  LONG_COMMAND("eval-dups-sbcds", &cmdEvalDupsSBCD, "Evaluate duplicates in spatial barcodes")

  LONG_COMMAND_GROUP("Tools for Analyzing Spatial Gene Expression (SGE)", NULL)
  LONG_COMMAND("make-sdge", &cmdMakeSpatialDGE, "Make spatial DGE files from STARsolo output")
  LONG_COMMAND("combine-sge", &cmdCombineSGE, "Combine multiple SGE files")
  LONG_COMMAND("dge2sdge", &cmdDGE2SDGE, "Convert DGE (from STARsolo) into SGE format (same as dge2sge)")
  LONG_COMMAND("dge2sge", &cmdDGE2SDGE, "Convert DGE (from STARsolo) into SGE format (same as dge2sdge)")
  LONG_COMMAND("sge2tsv", &cmdSGE2TSV, "Convert SGE (from sge2sdge) into plain TSV format (unsorted)")
  LONG_COMMAND("convert-sge", &cmdConvertSGE, "Convert SGE into plain various formats")
  LONG_COMMAND("segment-sge", &cmdSegmentSGE, "Segment SGE based on custom, grid, or hexagonal masks")
  LONG_COMMAND("subset-sge", &cmdSubsetSGE, "Subset Spatial SGE based on bounding box")
  LONG_COMMAND("combine-sge", &cmdCombineSGE, "Combine multiple SGE files")
  LONG_COMMAND("match-tag", &cmdMatchTag, "Match tags from FASTQ")
  LONG_COMMAND("search-tag", &cmdSearchTag, "Search tags from FASTQ")
  LONG_COMMAND("merge-matched-tags", &cmdMergeMatchedTags, "Merged matched tags generated by match-tag command")
  LONG_COMMAND("join-pixel-tsv", &cmdJoinPixelTSV, "Join pixel-level output from FICTURE with raw transcript-level TSV files")

  LONG_COMMAND_GROUP("Tools for FASTQ processing", NULL)
  LONG_COMMAND("filter-fastqs", &cmdFilterFASTQs, "Filter FASTQs based on a pattern")
  LONG_COMMAND("reformat-fastqs", &cmdReformatFASTQs, "Reformat FASTQs to be ready for STARsolo alignment")
  LONG_COMMAND("reformat-fastqs-mt", &cmdReformatFASTQsMT, "Reformat FASTQs to be ready for STARsolo alignment (multithreaded/not working)")
  LONG_COMMAND("filter-common-barcodes", &cmdFilterCommonBarcodes, "Filter FASTQ files based on the frequency of barcodes")
  LONG_COMMAND("stratify-fastq-by-barcodes", &cmdStratifyFASTQByBarcodes, "Stratify FASTQ files by spatial barcodes and UMIs")
  LONG_COMMAND("custom-demux-fastq", &cmdCustomDemuxFASTQ, "Demultiplex FASTQ files based in a customized manner")

  LONG_COMMAND_GROUP("Tools for Simple image processing", NULL)
  LONG_COMMAND("draw-xy", &cmdDrawXY, "Draw the image of points in 2D space")
  LONG_COMMAND("draw-3way", &cmdDraw3way, "Draw the 3-way image from the output of sttools pipeline")
  LONG_COMMAND("draw-sge", &cmdDrawSGE, "Draw the image of spatial gene expression (SGE) data")

  LONG_COMMAND_GROUP("Miscellaneous Tools", NULL)
  LONG_COMMAND("hist", &cmdHist, "Create a text-based histogram")

  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));

  if (argc < 2)
  {
    printf("[spatula] -- spatial transcriptomics utilities\n\n");
    fprintf(stderr, " Copyright (c) 2022-2024 by Hyun Min Kang\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 http://www.apache.org/licenses/\n\n");
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n", argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n", argv[0]);
    cl.Status();
    return 1;
  }
  else
  {
    if (strcmp(argv[1], "--help") == 0)
    {
      cl.HelpMessage();
    }
    else
    {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
