#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

/////////////////////////////////////////////////////////////////////////
// dge-barcode-match : Match a pair of barcode libraries
////////////////////////////////////////////////////////////////////////
int32_t cmdMakeSpatialDGE(int32_t argc, char** argv) {
  std::string bcdf;
  std::string ref_bcdf;
  std::string tgt_bcdf;
  std::string tgt_fastqf;
  std::string outf;
  int32_t len_bcd = 0;
  int32_t lunit_kmer = 0;
  int32_t nunit_kmer = 0;
  int32_t chunk_size = 100000;
  int32_t max_dist = -1;
  bool skip_no_match = false;
  bool skip_xtra = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("ref", &ref_bcdf, "Reference barcode dictionary")
    LONG_STRING_PARAM("tgt", &tgt_bcdf, "Target barcode to map to the reference")
    LONG_STRING_PARAM("fq", &tgt_fastqf, "Target FASTQ file to map to the reference, trimming and ignoring qualities")    

    LONG_PARAM_GROUP("Parameters for spaced k-mers", NULL)    
    LONG_INT_PARAM("len", &len_bcd, "Expected length of barcode")
    LONG_INT_PARAM("lunit", &lunit_kmer, "Expected length of kmer units (must be a multiple of 4)")
    LONG_INT_PARAM("nunit", &nunit_kmer, "Expected number of kmer units to make a single key to declare a match")
    LONG_INT_PARAM("chunk-size", &chunk_size, "Size of chunks to store spaced kmers")
    LONG_INT_PARAM("max-dist", &max_dist, "Maximum edit distance to allow")
    LONG_PARAM("skip-no-match", &skip_no_match, "Skip reads that does not match")
    LONG_PARAM("skip-xtra", &skip_xtra, "Skip writing extra fields even if it is available")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outf,"Output file name")
    LONG_INT_PARAM("verbose", &globalVerbosityThreshold, "Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( ref_bcdf.empty() || ( tgt_bcdf.empty() && tgt_fastqf.empty() ) || outf.empty() )
    error("Missing required option(s) : --ref, --tgt, --out");

  if ( len_bcd * lunit_kmer * nunit_kmer == 0 )
    error("Missing required options(s) : --len, --lunit, --nunit");

  if ( !tgt_bcdf.empty() && !tgt_fastqf.empty() )
    error("Cannot have both --tgt and --fq options");
  
  return 0;
}
