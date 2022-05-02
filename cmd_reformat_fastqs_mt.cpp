#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include <ctime>
#include <cstdio>
#include <set>
#include <unordered_set>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_LINE 65536
/////////////////////////////////////////////////////////////////////////
// reformat-fastq : Reformat FASTQ files before STARsolo alignment
////////////////////////////////////////////////////////////////////////
int32_t cmdReformatFASTQsMT(int32_t argc, char** argv) {
  std::string fq1f;
  std::string fq2f;
  std::string out1f;
  std::string out2f;  
  int32_t lenSBCD = 30;  // length of spatial barcode
  int32_t lenUMI = 9;  // length of UMI barcode in R2 (copied into R1)
  int32_t lenR2 = 101; // length of Read 2 after trimming.
  std::vector<std::string> matchTsvs;  // restrict to the matched spatial barcode
  int32_t lenMatch = 27; // length of perfect match requested

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("fq1", &fq1f, "FASTQ file for read 1 (plain file)")
    LONG_STRING_PARAM("fq2", &fq2f, "FASTQ file for read 2 (plain file)")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_MULTI_STRING_PARAM("match-tsv", &matchTsvs, ".matched.tsv.gz file(s) to filter FASTQ file based on spatial barcodes")
    LONG_INT_PARAM("len-match", &lenMatch, "Length of perfect match required with 2nd-seq FASTQ (27 or less)")    
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out1", &out1f, "Output FASTQ file for read1 (plain file)")
    LONG_STRING_PARAM("out2", &out2f, "Output FASTQ file for read2 (plain file)")
    LONG_INT_PARAM("len-sbcd", &lenSBCD, "Length of spatial barcode (Read 1)")
    LONG_INT_PARAM("len-umi", &lenUMI, "Length of UMI or randomer (Read 2)")
    LONG_INT_PARAM("len-r2", &lenR2, "Length of Read 2 to trim (including randomers)")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( fq1f.empty() || fq2f.empty() || out1f.empty() || out2f.empty() ) {
    error("Missing required options --fq1, --fq2, --out1, --out2");
  }

  if ( lenMatch > 27 )
    error("--len-match (%d) has to be 27 or less in the current implementation", lenMatch);

  // Read the pairs of FASTQ files
  FILE* rf1 = fopen(fq1f.c_str(), "r");
  if ( rf1 == NULL ) error("Cannot read %s", fq1f.c_str());
  
  FILE* rf2 = fopen(fq2f.c_str(), "r");
  if ( rf2 == NULL ) error("Cannot read %s", fq2f.c_str());

  notice("Reading input FASTQ file pairs");
  char name1[MAX_LINE];
  char name2[MAX_LINE];
  char seq1[MAX_LINE];
  char seq2[MAX_LINE];
  char dummy1[MAX_LINE];
  char dummy2[MAX_LINE];
  char qual1[MAX_LINE];
  char qual2[MAX_LINE];

  std::set<uint64_t> matchSet;
  uint64_t nmatch = 0;  
  for(int32_t i=0; i < (int32_t)matchTsvs.size(); ++i) {
    notice("Loading matched barcodes from %s", matchTsvs[i].c_str() );
    uint64_t sbcd;
    tsv_reader mtr(matchTsvs[i].c_str());
    while( mtr.read_line() > 0 ) {
      sbcd = seq2nt5(mtr.str_field_at(0),lenMatch);
      matchSet.insert(sbcd);
      ++nmatch;
    }
    notice("Loaded a total of %llu barcodes, %zu after removing redundant records", nmatch, matchSet.size());
    mtr.close();
  }

  FILE* wf1 = fopen(out1f.c_str(), "w");
  if ( wf1 == NULL ) error("Cannot write to %s", out1f.c_str());

  FILE* wf2 = fopen(out2f.c_str(), "w");  
  if ( wf2 == NULL ) error("Cannot write to %s", out2f.c_str());  
  
  notice("Writing to FASTQ file pairs %s and %s", out1f.c_str(), out2f.c_str());
  
  uint64_t nrecs = 0;
  bool skip = false;
  uint32_t nskip = 0;
  while( fgets(name1, MAX_LINE, rf1) != NULL ) {
    if ( ++nrecs % 10000000 == 0 ) 
      notice("Processing %llu records", nrecs);

    // parse the readname
    if ( fgets(name2, MAX_LINE, rf2) == NULL )
      error("Unexpected EOF in FASTQ file %s", fq2f.c_str());

    // read the sequence reads
    if ( fgets(seq1, MAX_LINE, rf1) == NULL )
      error("Unexpected EOF in FASTQ file %s", fq1f.c_str());
    if ( fgets(seq2, MAX_LINE, rf2) == NULL )
      error("Unexpected EOF in FASTQ file %s", fq2f.c_str());

    if ( !matchTsvs.empty() && matchSet.find(seq2nt5(seq1, lenMatch)) == matchSet.end() ) {
      skip = true;
      ++nskip;
    }
    else skip = false;

    if ( !skip ) {
      // write both read names without modifications
      fprintf(wf1, "%s", name1);
      fprintf(wf2, "%s", name2);

      // write modified sequence reads
      strncpy(seq1 + lenSBCD, seq2, lenUMI);
      seq1[lenSBCD+lenUMI] = '\0';
      seq2[lenR2] = '\0';
      int32_t lseq2 = strlen(seq2);
      if ( seq2[lseq2-1] == '\n' )
        seq2[lseq2-1] = '\0';
      fprintf(wf1,"%s\n", seq1);
      fprintf(wf2,"%s\n", seq2);
    }

    // read the dummy lines
    if ( fgets(dummy1, MAX_LINE, rf1) == NULL )
      error("Unexpected EOF in FASTQ file %s", fq1f.c_str());
    if ( fgets(dummy2, MAX_LINE, rf2) == NULL )
      error("Unexpected EOF in FASTQ file %s", fq2f.c_str());
    
    if ( dummy1[0] != '+' ) error("Unexpected line in FASTQ file %s", out1f.c_str());
    if ( dummy2[0] != '+' ) error("Unexpected line in FASTQ file %s", out2f.c_str());
    
    if ( !skip ) {
      fprintf(wf1, "%s", dummy1);
      fprintf(wf2, "%s", dummy2);
    }
    
    // read the quality strings
    if ( fgets(qual1, MAX_LINE, rf1) == NULL )
      error("Unexpected EOF in FASTQ file %s", out1f.c_str());
    if ( fgets(qual2, MAX_LINE, rf2) == NULL )
      error("Unexpected EOF in FASTQ file %s", out1f.c_str());      
    if ( !skip ) {
      // write modified sequence reads
      strncpy(qual1 + lenSBCD, qual2, lenUMI);
      qual1[lenSBCD+lenUMI] = '\0';
      qual2[lenR2] = '\0';
      int32_t lqual2 = strlen(qual2);
      if ( qual2[lqual2-1] == '\n' )
        qual2[lqual2-1] = '\0';      
      fprintf(wf1,"%s\n", qual1);
      fprintf(wf2,"%s\n", qual2);
    }
  }
  fclose(rf1);
  fclose(rf2);
  fclose(wf1);
  fclose(wf2);

  notice("Finished processing %llu records, skipping %llu.", nrecs-nskip, nskip);  

  notice("Analysis finished");
  
  return 0;
}
