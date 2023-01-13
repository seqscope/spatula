#include "spatula.h"
#include "qgenlib/commands.h"
#include "qgenlib/qgen_utils.h"

int32_t cmdMatchTag(int32_t argc, char** argv);
int32_t cmdSearchTag(int32_t argc, char** argv);
int32_t cmdMakeSpatialDGE(int32_t argc, char** argv);
int32_t cmdMatchSpatialBarcodes(int32_t argc, char** argv);
int32_t cmdMergeMatchedTags(int32_t argc, char** argv);
int32_t cmdBuildSpatialBarcodeDict(int32_t argc, char** argv);
int32_t cmdReformatFASTQs(int32_t argc, char** argv);
int32_t cmdFilterFASTQs(int32_t argc, char** argv);
int32_t cmdReformatFASTQsMT(int32_t argc, char** argv);
int32_t cmdDGE2SDGE(int32_t argc, char** argv);
 
int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND_GROUP("Spatial Transcriptomics Analysis Tools", NULL)
    LONG_COMMAND("make-sdge", &cmdMakeSpatialDGE, "Make spatial DGE files from STARsolo output")
    LONG_COMMAND("build-sbcds", &cmdBuildSpatialBarcodeDict, "Create spatial barcode dictionary based from HDMI FASTQ arrays")
    LONG_COMMAND("match-sbcds", &cmdMatchSpatialBarcodes, "Match the FASTQ file containing spatial barcodes with the spatial barcode dictionary")
    LONG_COMMAND("filter-fastqs", &cmdFilterFASTQs, "Filter FASTQs based on a pattern")    
    LONG_COMMAND("reformat-fastqs", &cmdReformatFASTQs, "Reformat FASTQs to be ready for STARsolo alignment")
    LONG_COMMAND("reformat-fastqs-mt", &cmdReformatFASTQsMT, "Reformat FASTQs to be ready for STARsolo alignment (multithreaded)")
    LONG_COMMAND("dge2sdge", &cmdDGE2SDGE, "Convert DGE (from STARsolo) into SDGE format")
    LONG_COMMAND("match-tag", &cmdMatchTag, "Match tags from FASTQ")
    LONG_COMMAND("search-tag", &cmdSearchTag, "Search tags from FASTQ")    
    LONG_COMMAND("merge-matched-tags", &cmdMergeMatchedTags, "Merged matched tags generated by match-tag command")        
  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));
  
  if ( argc < 2 ) {
    printf("[spatula] -- spatial transcriptomics utilities\n\n");
    fprintf(stderr, " Copyright (c) 2022 by Hyun Min Kang\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 http://www.apache.org/licenses/\n\n");    
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n",argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n",argv[0]);        
    cl.Status();
    return 1;
  }
  else {
    if ( strcmp(argv[1],"--help") == 0 ) {
      cl.HelpMessage();
    }
    else {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
