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
// filter-fastq : Filter FASTQ files based on certain patterns
////////////////////////////////////////////////////////////////////////
int32_t cmdFilterFASTQs(int32_t argc, char** argv) {
  std::string fq1f;
  std::string fq2f;
  std::string out1f;
  std::string out2f;
  std::vector<std::string> pat1s;
  std::vector<std::string> pat2s;
  int32_t min_mismatch = 0;
  bool remove = false;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("fq1", &fq1f, "FASTQ file for read 1")
    LONG_STRING_PARAM("fq2", &fq2f, "FASTQ file for read 2")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_MULTI_STRING_PARAM("pat1", &pat1s, "IUPAC pattern to match for Read 1 (1+ match required)")
    LONG_MULTI_STRING_PARAM("pat2", &pat2s, "IUPAC pattern to match for Read 2 (1+ match required)")
    LONG_INT_PARAM("min-mismatch", &min_mismatch, "Minimum number of mismatches allowed per pattern")
    LONG_PARAM("remove", &remove, "Remove the matching sequence from the FASTQ file instead of keeping it")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out1", &out1f, "Output FASTQ file for read1")
    LONG_STRING_PARAM("out2", &out2f, "Output FASTQ file for read2")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( fq1f.empty() || fq2f.empty() || out1f.empty() || out2f.empty() ) {
    error("Missing required options --fq1, --fq2, --out1, --out2");
  }

  // Read the pairs of FASTQ files
  htsFile* hf1 = hts_open(fq1f.c_str(), "r");
  if ( hf1 == NULL ) error("Cannot read %s", fq1f.c_str());
  
  htsFile* hf2 = hts_open(fq2f.c_str(), "r");
  if ( hf2 == NULL ) error("Cannot read %s", fq2f.c_str());

  notice("Reading FASTQ file pairs %s and %s", fq1f.c_str(), fq2f.c_str());  

  htsFile* wf1 = hts_open(out1f.c_str(), out1f.substr(out1f.size()-3).compare(".gz") == 0 ? "wz" : "w");
  if ( wf1 == NULL ) error("Cannot write to %s", out1f.c_str());

  htsFile* wf2 = hts_open(out2f.c_str(), out1f.substr(out2f.size()-3).compare(".gz") == 0 ? "wz" : "w");  
  if ( wf2 == NULL ) error("Cannot write to %s", out2f.c_str());

  
  notice("Writing to FASTQ file pairs %s and %s", out1f.c_str(), out2f.c_str());
  
  kstring_t str1; str1.l = str1.m = 0; str1.s = NULL;
  kstring_t str2; str2.l = str2.m = 0; str2.s = NULL;
  kstring_t name1; name1.l = name1.m = 0; name1.s = NULL;
  kstring_t name2; name2.l = name2.m = 0; name2.s = NULL;    
  uint64_t nrecs = 0;
  int32_t lstr1, lstr2, lname1, lname2;
  char buf1[65536];
  bool skip = false, pass1 = false, pass2 = false;
  uint32_t nskip = 0, npass = 0;

  std::vector<int32_t> lpat1s, lpat2s;
  for(int32_t i=0; i < (int32_t)pat1s.size(); ++i) {
    lpat1s.push_back((int32_t)pat1s[i].size());
  }
  for(int32_t i=0; i < (int32_t)pat2s.size(); ++i) {
    lpat2s.push_back((int32_t)pat2s[i].size());
  }
  while( (lname1 = hts_getline(hf1, KS_SEP_LINE, &name1)) > 0 ) {
    if ( ++nrecs % 1000000 == 0 )
      notice("Processing %llu records (%.5lf passed)", nrecs, (double)npass/(double)nrecs);

    // parse the readname
    lname2 = hts_getline(hf2, KS_SEP_LINE, &name2);
    if ( lname2 == 0 ) error("Unexpected EOF in FASTQ file %s", out2f.c_str());

    // read the sequence reads
    lstr1 = hts_getline(hf1, KS_SEP_LINE, &str1);
    if ( lstr1 == 0 ) error("Unexpected EOF in FASTQ file %s", out1f.c_str());    
    lstr2 = hts_getline(hf2, KS_SEP_LINE, &str2);
    if ( lstr2 == 0 ) error("Unexpected EOF in FASTQ file %s", out2f.c_str());

    pass1 = pat1s.empty();
    pass2 = pat2s.empty();

    for(int32_t i=0; i < (int32_t)pat1s.size(); ++i) {
      if ( seq_iupac_mismatch(str1.s, pat1s[i].c_str(), lpat1s[i]) <= min_mismatch ) {
        pass1 = true;
        break;
      }
    }

    for(int32_t i=0; i < (int32_t)pat2s.size(); ++i) {
      if ( seq_iupac_mismatch(str2.s, pat2s[i].c_str(), lpat2s[i]) <= min_mismatch ) {
        pass2 = true;
        break;
      }
    }

    skip = !(pass1 && pass2);
    if ( remove ) { skip = !skip; } // flip the filter if --remove was set
    if ( !skip ) {
      // write both read names without modifications
      hprintf(wf1, "%s\n", name1.s);
      hprintf(wf2, "%s\n", name2.s);
      hprintf(wf1, "%s\n", str1.s);
      hprintf(wf2, "%s\n", str2.s);
      ++npass;
    }
    else {
      ++nskip;
    }

    // read the dummy lines
    lstr1 = hts_getline(hf1, KS_SEP_LINE, &str1);
    if ( lstr1 == 0 ) error("Unexpected EOF in FASTQ file %s", fq1f.c_str());    
    lstr2 = hts_getline(hf2, KS_SEP_LINE, &str2);
    if ( lstr2 == 0 ) error("Unexpected EOF in FASTQ file %s", fq2f.c_str());
    if ( !skip ) {
      hprintf(wf1, "%s\n", str1.s);
      hprintf(wf2, "%s\n", str2.s);
    }    

    // read the quality strings lines
    lstr1 = hts_getline(hf1, KS_SEP_LINE, &str1);
    if ( lstr1 == 0 ) error("Unexpected EOF in FASTQ file %s", fq1f.c_str());    
    lstr2 = hts_getline(hf2, KS_SEP_LINE, &str2);
    if ( lstr2 == 0 ) error("Unexpected EOF in FASTQ file %s", fq2f.c_str());
    if ( !skip ) {
      hprintf(wf1, "%s\n", str1.s);
      hprintf(wf2, "%s\n", str2.s);
    }    
  }
  hts_close(hf1);
  hts_close(hf2);
  hts_close(wf1);
  hts_close(wf2);

  notice("Finished processing %llu records, skipping %llu.", npass, nskip);  

  free(str1.s);
  free(str2.s);
  free(name1.s);
  free(name2.s);    

  notice("Analysis finished");
  
  return 0;
}
