#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include <ctime>
#include <set>
#include <unordered_set>
#include <sys/stat.h>
#include <sys/types.h>

/////////////////////////////////////////////////////////////////////////
// reformat-fastq : Reformat FASTQ files before STARsolo alignment
////////////////////////////////////////////////////////////////////////
int32_t cmdReformatFASTQs(int32_t argc, char** argv) {
  std::string fq1f;
  std::string fq2f;
  std::string out1f;
  std::string out2f;  
  int32_t skipSBCD = 0;  // base to skip the first reads in Read 1
  int32_t lenSBCD = 30;  // length of spatial barcode
  int32_t lenUMI = 9;  // length of UMI barcode in R2 (copied into R1)
  int32_t lenR2 = 101; // length of Read 2 after trimming.
  std::string suffixR1(".R1.fasta.gz"); // suffix for read 1
  std::string suffixR2(".R2.fasta.gz"); // suffix for read 2
  std::vector<std::string> matchTsvs;  // restrict to the matched spatial barcode
  int32_t lenMatch = 27; // length of perfect match requested

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("fq1", &fq1f, "FASTQ file for read 1")
    LONG_STRING_PARAM("fq2", &fq2f, "FASTQ file for read 2")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_MULTI_STRING_PARAM("match-tsv", &matchTsvs, ".matched.tsv.gz file(s) to filter FASTQ file based on spatial barcodes")
    LONG_INT_PARAM("len-match", &lenMatch, "Length of perfect match required with 2nd-seq FASTQ (27 or less)")    
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out1", &out1f, "Output FASTQ file for read1")
    LONG_STRING_PARAM("out2", &out2f, "Output FASTQ file for read2")
    LONG_INT_PARAM("skip-sbcd", &skipSBCD, "Skip first bases of spatial barcode (Read 1)")
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
  htsFile* hf1 = hts_open(fq1f.c_str(), "r");
  if ( hf1 == NULL ) error("Cannot read %s", fq1f.c_str());
  
  htsFile* hf2 = hts_open(fq2f.c_str(), "r");
  if ( hf2 == NULL ) error("Cannot read %s", fq2f.c_str());

  notice("Reading FASTQ file pairs %s and %s", fq1f.c_str(), fq2f.c_str());  

  htsFile* wf1 = hts_open(out1f.c_str(), out1f.substr(out1f.size()-3).compare(".gz") == 0 ? "wz" : "w");
  if ( wf1 == NULL ) error("Cannot write to %s", out1f.c_str(), suffixR1.c_str());

  htsFile* wf2 = hts_open(out2f.c_str(), out1f.substr(out2f.size()-3).compare(".gz") == 0 ? "wz" : "w");  
  if ( wf2 == NULL ) error("Cannot write to %s", out2f.c_str(), suffixR2.c_str());

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
  
  notice("Writing to FASTQ file pairs %s and %s", out1f.c_str(), out2f.c_str());
  
  kstring_t str1; str1.l = str1.m = 0; str1.s = NULL;
  kstring_t str2; str2.l = str2.m = 0; str2.s = NULL;
  kstring_t name1; name1.l = name1.m = 0; name1.s = NULL;
  kstring_t name2; name2.l = name2.m = 0; name2.s = NULL;    
  uint64_t nrecs = 0;
  int32_t lstr1, lstr2, lname1, lname2;
  char buf1[65536];
  bool skip = false;
  uint32_t nskip = 0;
  while( (lname1 = hts_getline(hf1, KS_SEP_LINE, &name1)) > 0 ) {
    if ( ++nrecs % 10000000 == 0 ) 
      notice("Processing %llu records", nrecs);

    // parse the readname
    lname2 = hts_getline(hf2, KS_SEP_LINE, &name2);
    if ( lname2 == 0 ) error("Unexpected EOF in FASTQ file %s", out2f.c_str());

    // read the sequence reads
    lstr1 = hts_getline(hf1, KS_SEP_LINE, &str1);
    if ( lstr1 == 0 ) error("Unexpected EOF in FASTQ file %s", out1f.c_str());    
    lstr2 = hts_getline(hf2, KS_SEP_LINE, &str2);
    if ( lstr2 == 0 ) error("Unexpected EOF in FASTQ file %s", out2f.c_str());

    if ( !matchTsvs.empty() && matchSet.find(seq2nt5(str1.s + skipSBCD, lenMatch)) == matchSet.end() ) {
      skip = true;
      ++nskip;
    }
    else skip = false;

    if ( !skip ) {
      // write both read names without modifications
      hprintf(wf1, "%s\n", name1.s);
      hprintf(wf2, "%s\n", name2.s);

      // write modified sequence reads
      strncpy(buf1, str1.s + skipSBCD, lenSBCD);
      strncpy(buf1 + lenSBCD, str2.s, lenUMI);
      buf1[lenSBCD+lenUMI] = '\0';
      str2.s[lenR2] = '\0';
      hprintf(wf1,"%s\n", buf1);
      hprintf(wf2,"%s\n", str2.s);
    }

    // read the dummy lines
    lstr1 = hts_getline(hf1, KS_SEP_LINE, &str1);
    if ( lstr1 == 0 ) error("Unexpected EOF in FASTQ file %s", out1f.c_str());    
    lstr2 = hts_getline(hf2, KS_SEP_LINE, &str2);
    if ( lstr2 == 0 ) error("Unexpected EOF in FASTQ file %s", out2f.c_str());
    if ( str1.s[0] != '+' ) error("Unexpected line in FASTQ file %s", out1f.c_str());
    if ( str2.s[0] != '+' ) error("Unexpected line in FASTQ file %s", out2f.c_str());
    if ( !skip ) {
      hprintf(wf1, "%s\n", str1.s);
      hprintf(wf2, "%s\n", str2.s);
    }
    
    // read the quality strings
    lstr1 = hts_getline(hf1, KS_SEP_LINE, &str1);
    if ( lstr1 == 0 ) error("Unexpected EOF in FASTQ file %s", out1f.c_str());    
    lstr2 = hts_getline(hf2, KS_SEP_LINE, &str2);
    if ( lstr2 == 0 ) error("Unexpected EOF in FASTQ file %s", out2f.c_str());
    if ( !skip ) {
      strncpy(buf1, str1.s + skipSBCD, lenSBCD);
      strncpy(buf1 + lenSBCD, str2.s, lenUMI);
      buf1[lenSBCD+lenUMI] = '\0';
      str2.s[lenR2] = '\0';
      hprintf(wf1,"%s\n", buf1);
      hprintf(wf2,"%s\n", str2.s);
    }
  }
  hts_close(hf1);
  hts_close(hf2);
  hts_close(wf1);
  hts_close(wf2);

  notice("Finished processing %llu records, skipping %llu.", nrecs-nskip, nskip);  

  free(str1.s);
  free(str2.s);
  free(name1.s);
  free(name2.s);    

  notice("Analysis finished");
  
  return 0;
}
