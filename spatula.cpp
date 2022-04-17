#include "spatula.h"
#include "qgenlib/commands.h"
#include "qgenlib/qgen_utils.h"

int32_t cmdMakeSpatialDGE(int32_t argc, char** argv);
int32_t cmdMatchSpatialBarcodes(int32_t argc, char** argv);
int32_t cmdBuildSpatialBarcodeDict(int32_t argc, char** argv);
 
int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND_GROUP("Spatial Transcriptomics Analysis Tools", NULL)
    LONG_COMMAND("make-sdge", &cmdMakeSpatialDGE, "Make spatial DGE files from STARsolo output")
    LONG_COMMAND("build-sbcds", &cmdBuildSpatialBarcodeDict, "Create spatial barcode dictionary based from HDMI FASTQ arrays")
    LONG_COMMAND("match-sbcds", &cmdMatchSpatialBarcodes, "Match the FASTQ file containing spatial barcodes with the spatial barcode dictionary")          
  //LONG_COMMAND("match-sbcds", &cmdMatchSpatialBarcodes, "Compute summary statistics between FASTQ and spatial barcodes")  
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
