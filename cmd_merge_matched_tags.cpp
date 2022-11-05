#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "file_utils.h"
#include <ctime>
#include <set>
#include <cstdint>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

#define MAX_NT5_LEN 27
// #define UINT64_MAX ULLLONG_MAX
// #define INT32_MAX INT_MAX

// A class to sync barcode, tag, and UMI
class bcd_tag_umi_sync {
public:
  int32_t nbatches;  // number of batches to merge
  std::vector<uint64_t> nt5_bcds;  // max-27bp spatial barcode
  std::vector<int32_t>  int_tags;  // integer tag
  std::vector<uint64_t> nt5_umis;  // max-27bp UMIs

  // keeps track of smallest [spatial barcode, tag, UMIs]
  uint64_t min_bcd; 
  int32_t  min_tag;
  uint64_t min_umi;
  bool bcd_updated;
  bool tag_updated;
  bool umi_updated;
  std::vector<int32_t> imin_bcds; // keeps track of batches with smallest spatial barcode
  std::vector<int32_t> imin_tags; // keeps track of batches with smallest tag given spatial barcode
  std::vector<int32_t> imin_umis; // keeps track of batches with smallest UMI given tag and spatial barcode

  int32_t bcd_len, umi_len;  // length of raw barcode and UMIs
  int32_t nt5_bcd_len, nt5_umi_len; // length of barcodes and UMIs used to uint64_t encoding

  // initialize with # batches
  bcd_tag_umi_sync(int32_t _nbatches) : nbatches(_nbatches), min_bcd(UINT64_MAX-1), min_tag(INT32_MAX-1), min_umi(UINT64_MAX-1), bcd_updated(false), tag_updated(false), umi_updated(false), bcd_len(0), umi_len(0), nt5_bcd_len(0), nt5_umi_len(0) {
    nt5_bcds.resize(nbatches);
    int_tags.resize(nbatches);
    nt5_umis.resize(nbatches);
  }  

  // add barcode, tag, umi sequences from a specific batch
  // updates nt5_bcds, int_tags, and nt5_umis
  // return true if any entry was updated. return false if everything remains the same
  bool set_entry(int32_t ibatch, const char* bcd_seq, int32_t int_tag, const char* umi_seq) {
    if ( bcd_len == 0 ) { // first time setting up
      bcd_len = strlen(bcd_seq ? bcd_seq : "");
      umi_len = strlen(umi_seq ? umi_seq : "");
      nt5_bcd_len = bcd_len > MAX_NT5_LEN ? MAX_NT5_LEN : bcd_len;
      nt5_umi_len = umi_len > MAX_NT5_LEN ? MAX_NT5_LEN : umi_len;
    }
    uint64_t new_bcd = bcd_seq ? seq2nt5(bcd_seq, nt5_bcd_len) : UINT64_MAX;
    uint64_t new_umi = umi_seq ? seq2nt5(umi_seq, nt5_umi_len) : UINT64_MAX;
    bool ret = ((nt5_bcds[ibatch] != new_bcd) || (int_tags[ibatch] != int_tag) || (nt5_umis[ibatch] != new_umi));
    if ( ret ) { // if anything was changed, update.
      nt5_bcds[ibatch] = new_bcd; 
      int_tags[ibatch] = int_tag;
      nt5_umis[ibatch] = new_umi;
    }
    return ret;
  }

  // identify new minimum barcodes
  bool update_imin_bcds() {
    // find the minimum bcds first
    uint64_t old_bcd = min_bcd;
    min_bcd = nt5_bcds[0];
    imin_bcds.clear();
    imin_bcds.push_back(0);
    for(int32_t i=1; i < nbatches; ++i)  {
      if ( nt5_bcds[i] < min_bcd ) { // update min_bcd
        min_bcd = nt5_bcds[i];
        imin_bcds.clear();
        imin_bcds.push_back(i);
      }
      else if ( nt5_bcds[i] == min_bcd ) { // ties
        imin_bcds.push_back(i);
      }
    }
    bcd_updated = old_bcd != min_bcd;
    return bcd_updated;
  }

  // identify new minimum tag
  bool update_imin_tags() {
    int32_t old_tag = min_tag;
    // assuming imin_bcds were identified, identy minimum tags
    min_tag = int_tags[imin_bcds[0]];
    imin_tags.clear();
    imin_tags.push_back(imin_bcds[0]);
    for(int32_t i=1; i < (int32_t)imin_bcds.size(); ++i) {
      if ( int_tags[imin_bcds[i]] < min_tag ) { // update min_tag
        min_tag = int_tags[imin_bcds[i]];
        imin_tags.clear();
        imin_tags.push_back(imin_bcds[i]);
      }
      else if ( int_tags[imin_bcds[i]] < min_tag ) {
        imin_tags.push_back(imin_bcds[i]);        
      }
    }
    tag_updated = old_tag != min_tag;
    return tag_updated;
  }

  // identify new minimum umis
  bool update_imin_umis() {
    uint64_t old_umi = min_umi;
    // assuming imin_bcds were identified, identy minimum umis
    min_umi = nt5_umis[imin_tags[0]];
    imin_umis.clear();
    imin_umis.push_back(imin_tags[0]);
    for(int32_t i=1; i < (int32_t)imin_tags.size(); ++i) {
      if ( nt5_umis[imin_tags[i]] < min_umi ) { // update min_tag
        min_umi = nt5_umis[imin_tags[i]];
        imin_umis.clear();
        imin_umis.push_back(imin_tags[i]);
      }
      else if ( nt5_umis[imin_tags[i]] < min_umi ) {
        imin_umis.push_back(imin_tags[i]);  
      }
    }
    umi_updated = old_umi != min_umi;
    return umi_updated;    
  }

  // update minimums, return false when reaching EOF
  bool update_imins() {
    update_imin_bcds();
    update_imin_tags();
    update_imin_umis();
    return min_bcd != UINT64_MAX;
  }
};

int32_t cmdMergeMatchedTags(int32_t argc, char** argv) {
  std::string listf;      // file containing the list of batches to merge
  std::string tagf;       // dictionary of tags
  std::string outdir;     // output directory

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input files", NULL)
    LONG_STRING_PARAM("list", &listf, "List of output files from match-tag")
    LONG_STRING_PARAM("tag", &tagf, "Dictionary file containing tag sequences")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outdir,"Output directory to store DGE file")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( listf.empty() || tagf.empty() || outdir.empty() ) {
    error("Missing required options --list, --tag, or --out");
  }

  // create the output directory first
  if ( makePath(outdir) ) {
    notice("Successfully created the directory %s", outdir.c_str());
  } else {
    notice("The directory %s already exists", outdir.c_str());    
  }

  // process the file list
  tsv_reader tr_list(listf.c_str());
  std::vector<std::string> filenames;
  while( tr_list.read_line() ) {
    filenames.push_back(tr_list.str_field_at(0));
  }
  tr_list.close();

  // process the tags
  tsv_reader tr_tags(tagf.c_str());
  std::vector<std::string> tag_ids;
  std::vector<std::string> tag_names;
  while( tr_tags.read_line() ) {
    tag_ids.push_back(tr_tags.str_field_at(0));
    tag_names.push_back(tr_tags.str_field_at(1));    
  }
  tr_tags.close();

  // open each batch and peek the first line
  std::vector<tsv_reader*> batch_trs;
  int32_t nbatches = (int32_t)filenames.size();

  bcd_tag_umi_sync btus(nbatches);
  for(int32_t i=0; i < nbatches; ++i) {
    batch_trs.push_back(new tsv_reader(filenames[i].c_str()));
    batch_trs[i]->read_line(); // peek one line for each
    btus.set_entry(i, batch_trs[i]->str_field_at(0), batch_trs[i]->int_field_at(1), batch_trs[i]->str_field_at(2)); 
  }
  htsFile* wbcd = hts_open((outdir + "/barcodes.tsv.gz").c_str(), "wz");
  htsFile* wftr = hts_open((outdir + "/features.tsv.gz").c_str(), "wz");
  htsFile* wread = hts_open((outdir + "/.tmp.reads.mtx").c_str(), "w");
  htsFile* wumi = hts_open((outdir + "/.tmp.umis.mtx").c_str(), "w");  

  for(int32_t i=0; i < (int32_t)tag_ids.size(); ++i) {
    hprintf(wftr, "%s\t%s\tAntibody_Tag\n", tag_ids[i].c_str(), tag_names[i].c_str());
  }
  hts_close(wftr);

  int32_t nreads = 0, numis = 0;
  bool updated = false;
  uint64_t ibcd = 0, ilines = 0;
  char bcd_seq[255], umi_seq[255];
  while ( btus.update_imins() ) { // keep reading as long as EOF has not reached yet
    // keep reading ties
    for(int32_t i=0; i < (int32_t)btus.imin_umis.size(); ++i) {
      int32_t j = btus.imin_umis[i];
      if (i == 0) { // make a copy of spatial barcode and UMI, which are not stored in bcd_tag_umi_sync
        strcpy(bcd_seq, batch_trs[j]->str_field_at(0));
        strcpy(umi_seq, batch_trs[j]->str_field_at(2));
        if ( btus.bcd_updated ) {  // if barcode is updated, write it out
          hprintf(wbcd, "%s\n", batch_trs[j]->str_field_at(0));
          ++ibcd;

          if ( ibcd % 1000000 == 0 )
            notice("Writing %llu barcodes...", ibcd);
        }
      }      
      updated = false;
      while ( !updated ) {
        ++nreads;
        if ( batch_trs[j]->read_line() ) {
          updated = btus.set_entry(j, batch_trs[j]->str_field_at(0), batch_trs[j]->int_field_at(1), batch_trs[j]->str_field_at(2));
        }
        else {
          updated = btus.set_entry(j, NULL, INT32_MAX, NULL);
        }
      }
    }
    // if either barcode or tag was updated, there are no UMI ties, and should be written
    if (!btus.bcd_updated && !btus.tag_updated) { // only UMI were updated. Do not write yet
      ++numis;
    }
    else {
      hprintf(wread, "%d %lld %d\n", btus.min_tag + 1, ibcd, nreads);
      hprintf(wumi, "%d %lld %d\n", btus.min_tag + 1, ibcd, numis+1);      
      ++ilines;
      nreads = numis = 0;
    }
  }
  hts_close(wread);
  hts_close(wumi);  
  hts_close(wbcd);

  htsFile* whdr = hts_open((outdir + "/.tmp.matrix.hdr").c_str(), "w");
  hprintf(whdr,"%%MatrixMarket matrix coordinate integer general\n%\n");
  hprintf(whdr, "%zu %llu %llu\n", tag_ids.size(), ibcd, ilines);
  hts_close(whdr);

  notice("Generating the merged matrix.mtx.gz file");  
  std::string cmd;
  catprintf(cmd, "cat %s %s | gzip -c > %s", (outdir + "/.tmp.matrix.hdr").c_str(), (outdir + "/.tmp.reads.mtx").c_str(), (outdir + "/reads.mtx.gz").c_str());  
  int32_t ret = system(cmd.c_str()); // run the commnad
  if ( (ret == -1) || (WEXITSTATUS(ret) == 127) ) {
    error("Error in running %s", cmd.c_str());
  }
  cmd.clear();
  catprintf(cmd, "cat %s %s | gzip -c > %s", (outdir + "/.tmp.matrix.hdr").c_str(), (outdir + "/.tmp.umis.mtx").c_str(), (outdir + "/umis.mtx.gz").c_str());  
  ret = system(cmd.c_str()); // run the commnad
  if ( (ret == -1) || (WEXITSTATUS(ret) == 127) ) {
    error("Error in running %s", cmd.c_str());
  }  
  if ( remove((outdir + "/.tmp.matrix.hdr").c_str()) != 0 )
    error("Cannot remove %s", (outdir + "/.tmp.matrix.hdr").c_str());
  if ( remove((outdir + "/.tmp.reads.mtx").c_str()) != 0 )
    error("Cannot remove %s", (outdir + "/.tmp.reads.mtx").c_str());
  if ( remove((outdir + "/.tmp.umis.mtx").c_str()) != 0 )
    error("Cannot remove %s", (outdir + "/.tmp.umis.mtx").c_str());  

  notice("Analysis finished");
  
  return 0;
}
