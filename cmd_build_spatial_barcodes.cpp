#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
extern "C" {
#include "htslib/kstring.h"
}
#include <vector>
#include <string>
#include <cstring>
#include <climits>
#include <map>

class lane_tile_t {
public:
  htsFile* hp;
  std::string filename;
  int32_t nrecords;
  int32_t nmatches;
  int32_t xmin;
  int32_t xmax;
  int32_t ymin;
  int32_t ymax;

  lane_tile_t(htsFile* _hp, const char* _fn)
    : hp(_hp), filename(_fn), nrecords(0), nmatches(0),
      xmin(INT_MAX), xmax(0), ymin(INT_MAX), ymax(0) {
    //notice("cons %x %s", _hp, _fn);
  }
};

void parse_illumina_readname(kstring_t* str, int32_t* lane, int32_t* tile, int32_t* x, int32_t* y) {
  int32_t* fields = NULL;
  int32_t nfields;
  fields = ksplit(str, ':', &nfields);
  if ( nfields < 7 )
    error("Cannot parse Readname '%s'. Observed only %d fields when separating by colon", str->s, nfields);
  *lane = atoi(str->s + fields[3]);
  *tile = atoi(str->s + fields[4]);
  *x = atoi(str->s + fields[5]);
  *y = atoi(str->s + fields[6]);
  free(fields);
}

void parse_salus_readname(kstring_t* str, int32_t* tile, int32_t* x, int32_t* y) {
  int32_t* fields = NULL;
  int32_t nfields;
  fields = ksplit(str, '_', &nfields);

  if ( nfields < 6 ) {
    error("Cannot parse Readname '%s'. Observed only %d fields when separating by underscore", str->s, nfields);
  }

  const char* tilestr = str->s + fields[nfields-6];
  if ( ( tilestr[0] == 'R' ) && ( tilestr[3] == 'C' ) ) { // conforms to the expected format
    *tile = 1000000 + (atoi(tilestr+1) * 1000) + atoi(tilestr+4);
  }
  else {
    error("Cannot parse Readname '%s'. Cannot find the tile information", str->s);
  }

  *x = (int)(atof(str->s + fields[nfields-2])*1000);
  *y = (int)(atof(str->s + fields[nfields-1])*1000);
  free(fields);
}

void parse_salus_global_readname(kstring_t* str, int32_t* tile, int32_t* x, int32_t* y) {
  int32_t* fields = NULL;
  int32_t nfields;
  fields = ksplit(str, '_', &nfields);

  if ( nfields < 6 ) {
    error("Cannot parse Readname '%s'. Observed only %d fields when separating by underscore", str->s, nfields);
  }

  *tile = 1;
  *x = (int)(atof(str->s + fields[nfields-2])*1000);
  *y = (int)(atof(str->s + fields[nfields-1])*1000);
  free(fields);
}

/////////////////////////////////////////////////////////////////////////
// dge-barcode-match : Match a pair of barcode libraries
////////////////////////////////////////////////////////////////////////
int32_t cmdBuildSpatialBarcodeDict(int32_t argc, char** argv) {
  std::string fastqf;
  std::string format;
  std::string platform("Illumina");
  std::string outf;
  int32_t force_lane = 0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("fq", &fastqf, "FASTQ file that contains spatial barcode")
    LONG_STRING_PARAM("format", &format, "Format of the HDMI array (DraI32, DMix32, DraI20, DPAGE32, HDMI20, HDMI30, DraI31)")
    LONG_STRING_PARAM("platform", &platform, "Platform of the sequencing data to determine the rule to parse the readnames (accepting Illumina, Salus, SalusGlobal)")
    LONG_INT_PARAM("force-lane", &force_lane, "Force a lane number. Required a positive value for Salus/SalusGlobal platforms")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outf,"Output directory name")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( fastqf.empty() || outf.empty() ) {
    error("Missing required options --fq, --out");
  }  

  int32_t hdmi_length;
  std::vector<std::string> hdmi_patterns;
  bool rev_comp = true;

  if ( format.compare("DraI32") == 0 ) {
    hdmi_length = 32;
    hdmi_patterns.push_back("NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN");
  }
  else if ( format.compare("DraI31") == 0 ) {
    hdmi_length = 31;
    hdmi_patterns.push_back("NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBN");
  }
  else if ( format.compare("T7-30") == 0 ) {
    hdmi_length = 30;
    hdmi_patterns.push_back("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
    rev_comp = false;
  }    
  else if ( format.compare("HDMI30") == 0 ) {
    hdmi_length = 30;
    hdmi_patterns.push_back("NNNNNBNNBNNBNNBNNBNNBNNBNNBNNB");
  }    
  else if ( format.compare("DraI20") == 0 ) {
    hdmi_length = 20;
    hdmi_patterns.push_back("NNNNNBNNBNNBNNBNNBNN");
  }
  else if ( format.compare("HDMI20") == 0 ) {
    hdmi_length = 20;
    hdmi_patterns.push_back("NNHNVNNVNVNVNVNVNVNN");
  }
  else if ( format.compare("DMix32") == 0 ) {
    hdmi_length = 32;
    hdmi_patterns.push_back("NNNNNBNNBNNBNNBNNBNNBNNBNNBNNSNN");
    hdmi_patterns.push_back("NNNBNNBNNBNNBNNBNNBNNBNNBNNBNNBN");
    hdmi_patterns.push_back("NNNNBNNBNNBNNBNNBNNBNNBNNBNNBNNB");
  }
  else if ( format.compare("XGMix32") == 0 ) {
    hdmi_length = 32;
    hdmi_patterns.push_back("NNHHNNHHNNHHNNHHNNHHNNHHNNHHNNHH");
    hdmi_patterns.push_back("DDNNDDNNDDNNDDNNDDNNDDNNDDNNDDNN");
  }
  else if ( format.compare("DPAGE32") == 0 ) {
    hdmi_length = 32;
    hdmi_patterns.push_back("NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN");
    hdmi_patterns.push_back("NNNNBNNBNNBNNBNNBNNBNNBNNBNNBNNB");
  }
  else if ( format.compare("XG3") == 0 ) {
    hdmi_length = 32;
    hdmi_patterns.push_back("DDNNDDNNDDNNDDNNDDNNDDNNDDNNDDNN");
    hdmi_patterns.push_back("NNHHNNHHNNHHNNHHNNHHNNHHNNHHNNHH");
  }
  else {
    error("Cannot recognize the HDMI format %s. Acceptable values are DraI32, DMix32, DraI20, HDMI20", format.c_str());
  }

  if ( hdmi_length != hdmi_patterns[0].size() ) {
    error("HDMI length %d does not match with %s", hdmi_length, hdmi_patterns[0].c_str());
  }

  htsFile* hp = hts_open(fastqf.c_str(), "r");

  kstring_t str; str.l = str.m = 0; str.s = NULL;
  uint64_t nrecs = 0;
  int32_t* fields = NULL;
  int32_t i, j, lstr, nfields, lane, tile, x, y, lseq, ldummy, lqual, mismatches, tmp;
  char tmpbuf[65535];
  std::map<std::string,lane_tile_t*> ltdict;
  std::map<std::string,lane_tile_t*>::iterator ltitr;
  lane_tile_t* ltval;
  
  if ( outf[outf.size()-1] != '/' ) 
    outf += "/";

  enum platform_t { ILLUMINA, SALUS, SALUS_GLOBAL };
  platform_t platform_type = ILLUMINA;
  if ( platform.compare("Illumina") == 0 ) {
    platform_type = ILLUMINA;
  }
  else if ( platform.compare("Salus") == 0 ) {
    platform_type = SALUS;
    if ( force_lane <= 0 )
      error("For Salus platform, you must specify a lane number using --force-lane option");
  }
  else if ( platform.compare("SalusGlobal") == 0 ) {
    platform_type = SALUS_GLOBAL;
    if ( force_lane <= 0 )
      error("For Salus platform, you must specify a lane number using --force-lane option");
  }
  else {
    error("Cannot recognize the platform %s. Acceptable values are Illumina, Salus, SalusGlobal", platform.c_str());
  }

  while( (lstr = hts_getline(hp, KS_SEP_LINE, &str)) > 0 ) { 
    if ( ++nrecs % 1000000 == 0 )
      notice("Reading %llu records from %s, %zu tiles observed so far", nrecs, fastqf.c_str(), ltdict.size());

    // separate the Readname by ':'
    switch(platform_type) {
    case ILLUMINA:
      parse_illumina_readname(&str, &lane, &tile, &x, &y);
      break;
    case SALUS:
      lane = force_lane;
      parse_salus_readname(&str, &tile, &x, &y);
      break;
    case SALUS_GLOBAL:
      lane = force_lane;
      parse_salus_global_readname(&str, &tile, &x, &y);
      break;
    };
    fields = ksplit(&str, ':', &nfields);
    if ( nfields < 7 )
      error("Cannot parse Readname '%s' in FASTQ file %s at record=%llu. Observed only %d fields when separating by colon", str.s, fastqf.c_str(), nrecs, nfields);

    lane = atoi(str.s + fields[3]);
    tile = atoi(str.s + fields[4]);
    x = atoi(str.s + fields[5]);
    y = atoi(str.s + fields[6]);

    // tmpbuf is the key
    snprintf(tmpbuf, 65535, "%d_%d", lane, tile);
    ltitr = ltdict.find(tmpbuf);

    if ( ltitr == ltdict.end() ) { // need to create a new filehandle
      std::string fn = outf + tmpbuf + ".sbcds.unsorted.tsv";
      htsFile* wf = hts_open(fn.c_str(), "w");
      if ( wf == NULL )
        error("Cannot open a file %s for writing", fn.c_str());
      ltval = new lane_tile_t(wf, fn.c_str());
      ltdict[tmpbuf] = ltval;
    }
    else {
      ltval = ltitr->second;
    }

    //notice("bar");        

    // read the sequence reads for line 4N+1
    lseq = hts_getline(hp, KS_SEP_LINE, &str);
    if ( lseq < hdmi_length )
      error("Cannot parse sequence reads in FASTQ file %s at record=%llu. Read length is too short (%d)", fastqf.c_str(), nrecs, lseq);
    str.s[hdmi_length] = '\0'; // trim the sequence file

    // calculate mismatches from the expected pattern
    mismatches = seq_iupac_mismatch(str.s, hdmi_patterns[0].c_str(), hdmi_length);
    for(i=1; i < (int32_t)hdmi_patterns.size(); ++i) {
      //error("foo");
      tmp = seq_iupac_mismatch(str.s, hdmi_patterns[i].c_str(), hdmi_length);
      if ( tmp < mismatches )
        mismatches = tmp;
    }

    // update values associated with each tile
    ltval->nrecords += 1;
    if ( mismatches == 0 )
      ltval->nmatches += 1;
    if ( ltval->xmin > x ) ltval->xmin = x;
    if ( ltval->xmax < x ) ltval->xmax = x;
    if ( ltval->ymin > y ) ltval->ymin = y;
    if ( ltval->ymax < y ) ltval->ymax = y;        
    
    if ( rev_comp )
      seq_revcomp(str.s, hdmi_length); // take reverse complement of the sequence
    hprintf(ltval->hp, "%s\t%d\t%d\t%d\t%d\t%d\n", str.s, lane, tile, x, y, mismatches);
    //notice("%s\t%d\t%d\t%d\t%d\t%d", str.s, lane, tile, x, y, mismatches);
    free(fields);    

    ldummy = hts_getline(hp, KS_SEP_LINE, &str);
    lqual  = hts_getline(hp, KS_SEP_LINE, &str);    

    //if ( nrecs % 10 == 0 ) break;  // for debugging
  }

  for(ltitr = ltdict.begin(); ltitr != ltdict.end(); ++ltitr) {
    hts_close(ltitr->second->hp);
  }

  notice("Finished reading %llu records from %s, %zu tiles observed in total", nrecs, fastqf.c_str(), ltdict.size());

  // write the manifest file
  htsFile* wf = hts_open( (outf + "manifest.tsv").c_str(), "w" );
  if ( wf == NULL )
      error("Cannot open a file %smanifest.tsv for writing", outf.c_str());
  hprintf(wf, "id\tfilepath\tbarcodes\tmatches\tmismatches\txmin\txmax\tymin\tymax\n");
  for(ltitr = ltdict.begin(); ltitr != ltdict.end(); ++ltitr) {
    lane_tile_t* p = ltitr->second;
    const char* key = ltitr->first.c_str();
    // we write as .sbcds.sorted.tsv.gz because we expect that the pipeline will sort the tsv file later on
    hprintf(wf, "%s\t%s.sbcds.sorted.tsv.gz\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", key, key, p->nrecords, p->nmatches, p->nrecords - p->nmatches, p->xmin, p->xmax, p->ymin, p->ymax); 
  }
  hts_close(wf);

  notice("Analysis Finished");
  
  return 0;
}
