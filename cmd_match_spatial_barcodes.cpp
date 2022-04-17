#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

// open all tiles
void open_tiles(dataframe_t& df, std::vector<std::string>& tiles, std::vector<tsv_reader*>& bcdfs) {
  // read manifest files
  if ( bcdfs.size() > 0 ) {
    for(int32_t i=0; i < (int32_t) bcdfs.size(); ++i) {
      delete bcdfs[i];
    }
    bcdfs.clear();
  }
  
  tiles = df.get_column("id");
  int32_t icol = df.get_colidx("fullpath");
  for(int32_t i=0; i < df.nrows; ++i) {
    bcdfs.push_back(new tsv_reader(df.get_str_elem(i, icol).c_str()));
    if ( bcdfs.back()->read_line() == 0 ) {
      error("ERROR: Observed an empty barcode file %s", df.get_str_elem(i,icol));      
    }    
  }

  /*
  tsv_reader idxf((bcddir + "manifest.tsv").c_str());
  if ( idxf.read_line() == 0 ) // read the header
    error("Cannot read the header of %smanifest.tsv", idxf);

  while( idxf.read_line() > 0 ) { // read each line
    std::string fn = bcddir + idxf.str_field_at(1);
    tiles.push_back(idxf.str_field_at(0));
    bcdfs.push_back(new tsv_reader(fn.c_str()));
    if ( bcdfs.back()->read_line() == 0 ) {
      error("ERROR: Observed an empty barcode file %s", fn.c_str());      
    }
  }
  idxf.close();*/
}

uint64_t count_matches(std::vector<uint64_t>& bseqs, dataframe_t& df, std::vector<uint64_t>& counts, int32_t match_len, htsFile* wmatch) {
//uint64_t count_matches(std::vector<std::string>& bseqs, std::string& bcddir, std::vector<uint64_t>& counts) {
  std::vector<std::string> tiles;
  std::vector<tsv_reader*> bcdfs;

  open_tiles(df, tiles, bcdfs);

  int32_t ntiles = (int32_t)tiles.size();
  if ( counts.empty() ) 
    counts.resize(ntiles, 0);
  int32_t len = strlen(bcdfs[0]->str_field_at(0));
  //if ( len != (int32_t) bseqs[0].size() )
  //  error("HDMI length %d does not match to the parameters %zu", len, bseqs[0].size());  
  if ( len < match_len )
    error("HDMI length %d does not match to the parameters %d", len, match_len);

  /*
  char zzz[65535];
  for(int32_t i=0; i < len; ++i) {
    zzz[i] = 'Z';
  }
  zzz[len] = '\0';
  */

  std::vector<uint64_t> tseqs(ntiles);
  //std::vector<std::string> tseqs;
  for(int32_t i=0; i < ntiles; ++i) {
    tseqs[i] = seq2nt5(bcdfs[i]->str_field_at(0),match_len);    
    //tseqs[i] = seq2bits(bcdfs[i]->str_field_at(0),len);
    //tseqs.push_back(bcdfs[i]->str_field_at(0));
    //notice("%d\t%032llu", i, tseqs[i]);
  }

  // sort the batch of sequences
  uint64_t batch_size = (uint64_t)bseqs.size();
  notice("Started sorting of %llu records", batch_size);
  std::sort(bseqs.begin(), bseqs.end());
  notice("Finished sorting of %llu records", batch_size);  
  uint64_t nseqs = (uint64_t)bseqs.size();

  uint64_t nmiss = 0;
  uint64_t ndups = 0;
  bool has_match;
  int32_t cmp;
  // count the sequences that matches
  for(uint64_t i=0; i < nseqs; ++i) {
    if (i % (batch_size / 20) == 0)
      notice("Processing %d records, nmiss = %llu, ndups = %llu, bseqs[i]=%032llu, tseqs[0]=%032llu", i, nmiss, ndups, bseqs[i], tseqs[0]);
    if ( ( i > 0 ) && ( bseqs[i] == bseqs[i-1] ) ) { // avoid double-counting duplicate spatial barcodes
      ++ndups;
      continue;
    }
    has_match = false;
    //const std::string& s = bseqs[i];
    uint64_t s = bseqs[i];
    for(int32_t j=0; j < ntiles; ++j) {
      cmp = s < tseqs[j] ? -1 : ( s == tseqs[j] ? 0 : 1 );
      //cmp = s.compare(tseqs[j]);
      while( cmp > 0 ) {
        //if ( bcdfs[j]->read_line() == 0 ) tseqs[j].assign(zzz);
        //else tseqs[j].assign(bcdfs[j]->str_field_at(0));
        if ( bcdfs[j]->read_line() == 0 ) tseqs[j] = UINT64_MAX;
        else tseqs[j] = seq2nt5(bcdfs[j]->str_field_at(0),match_len);
        //else tseqs[j] = seq2bits(bcdfs[j]->str_field_at(0),len);        
        //cmp = s.compare(tseqs[j]);
        cmp = s < tseqs[j] ? -1 : ( s == tseqs[j] ? 0 : 1 );        
      }
      if ( cmp == 0 ) {
        has_match = true;
        hprintf(wmatch,"%s", bcdfs[j]->str_field_at(0));
        for(int32_t k=1; k < bcdfs[j]->nfields; ++k) 
          hprintf(wmatch,"\t%s", bcdfs[j]->str_field_at(k));
        hprintf(wmatch,"\n");
        ++counts[j];
      }
    }
    if ( !has_match ) {
      ++nmiss;
    }
  }
  
  for(int32_t i=0; i < ntiles; ++i) {
    delete bcdfs[i];
  }

  return nmiss;
}

/////////////////////////////////////////////////////////////////////////
// match-sbceds : Match spatial barcodes
////////////////////////////////////////////////////////////////////////
int32_t cmdMatchSpatialBarcodes(int32_t argc, char** argv) {
  std::string fastqf;
  std::string bcddir;
  std::string outprefix;
  int32_t batch_size = 300000000; // number of records to be searched for at once
  //int32_t hdmi_len = 32;
  int32_t match_len = 27;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("fastq", &fastqf, "FASTQ file read 1 containing 2nd-seq spatial barcode")
    LONG_STRING_PARAM("bcd", &bcddir, "Spatial barcode dictionary generated from 'build-sbcds' command")
    LONG_INT_PARAM("batch", &batch_size, "Size of a single batch")
    LONG_INT_PARAM("match-len", &match_len, "Length of HDMI spatial barcodes to require perfect matches")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outprefix,"Output prefix (index.tsv, matches.tsv.gz)")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( fastqf.empty() || bcddir.empty() || outprefix.empty() ) {
    error("Missing required options --fq, --bcd, --out");
  }

  // make spatial barcode directory end with '/'
  if ( bcddir[bcddir.size()-1] != '/' ) 
    bcddir += "/";

  // read the manifest file
  dataframe_t df((bcddir + "manifest.tsv").c_str());
  // add a column to indicate the full path
  int32_t icol = df.add_empty_column("fullpath");
  int32_t jcol = df.get_colidx("filepath");
  for(int32_t i=0; i < df.nrows; ++i) {
    df.set_str_elem((bcddir + df.get_str_elem(i, jcol)).c_str(), i, icol);
  }

  htsFile* wmatch = hts_open((outprefix + ".matched.tsv.gz").c_str(), "wz");
  if ( wmatch == NULL )
    error("Cannot open %s.matches.tsv.gz for writing", outprefix.c_str());

  // Read the 2nd-seq FASTQ
  htsFile* hp = hts_open(fastqf.c_str(), "r");
  notice("Reading FASTQ file %s", fastqf.c_str());  
  std::vector<uint64_t> bseqs;      // store the 2nd-sequence into an array after conveting to 64bit integer
  //std::vector<std::string> bseqs; // an exact thing to do is to store actual strings
  std::vector<uint64_t> counts;
  kstring_t str; str.l = str.m = 0; str.s = NULL;
  uint64_t nrecs = 0, ibatch = 0;
  uint64_t nmiss = 0;
  int32_t lstr, lseq, ldummy, lqual;
  while( (lstr = hts_getline(hp, KS_SEP_LINE, &str)) > 0 ) {
    if ( nrecs % 10000000 == 0 ) 
      notice("Processing %d records from the FASTQ file %s", nrecs, fastqf.c_str());      
    
    // read the sequence reads for line 4N+1
    lseq = hts_getline(hp, KS_SEP_LINE, &str);
    if ( lseq < match_len )
      error("Cannot parse Readname in FASTQ file %s at record=%llu. Read length is too short (%d)", fastqf.c_str(), nrecs, lseq);
    
    bseqs.push_back(seq2nt5(str.s,match_len));
    //bseqs.push_back(seq2bits(str.s,match_len));    
    //bseqs.push_back(str.s);
    ldummy = hts_getline(hp, KS_SEP_LINE, &str);
    lqual  = hts_getline(hp, KS_SEP_LINE, &str);

    if ( ++nrecs % batch_size == 0 ) {
      notice("Processing batch %d of %d sequences", ++ibatch, batch_size);
      //nmiss += count_matches(bseqs, bcddir, counts);
      nmiss += count_matches(bseqs, df, counts, match_len, wmatch);
      bseqs.clear();
    }
  }
  if ( bseqs.size() > 0 ) {
    notice("Processing the last batch %d containing %zu sequences", ++ibatch, bseqs.size());
    //nmiss += count_matches(bseqs, bcddir, counts);
    nmiss += count_matches(bseqs, df, counts, match_len, wmatch);
    bseqs.clear();
  }
  hts_close(wmatch);  

  // read manifest files
  htsFile* wf = hts_open((outprefix + ".counts.tsv").c_str(), "w");
  hprintf(wf, "id\tfilepath\tbarcodes\tcounts\n");
  for( int32_t i=0; i < df.nrows; ++i) {
    hprintf(wf, "%s", df.get_str_elem(i, "id").c_str());
    hprintf(wf, "\t%s", df.get_str_elem(i, "filepath").c_str());
    hprintf(wf, "\t%llu", df.get_uint64_elem(i, "barcodes"));
    hprintf(wf, "\t%llu\n", counts[i]);
  }
  hts_close(wf);
  
  notice("Total = %llu, Missed = %llu (%.5f)", nrecs, nmiss, (double)nmiss/(double)nrecs);  

  notice("Analysis finished");
  
  return 0;
}
