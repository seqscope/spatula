#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

#define MAX_NT5_LEN 27

int32_t cmdSearchTag(int32_t argc, char** argv) {
  std::string fq1f;
  std::string fq2f;
  std::string tagf;
  std::string smatchf;
  std::string bcddir;
  std::string outprefix;
  std::string umipos;
  std::string bcdpos = "1-32";
  //std::string tagpos = "1-15";
  int32_t batch_size = 10000000;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input files", NULL)
    LONG_STRING_PARAM("fq1", &fq1f, "FASTQ file read 1 containing 2nd-seq spatial barcode")
    LONG_STRING_PARAM("fq2", &fq2f, "FASTQ file read 2 containing 2nd-seq spatial barcode")
    LONG_STRING_PARAM("tag", &tagf, "Dictionary file containing tag sequences to search for")
    LONG_STRING_PARAM("smatch", &smatchf, "Output from smatch step that contains spatial barcodes to be used for whitelist")

    LONG_PARAM_GROUP("Positions in reads", NULL)
    LONG_STRING_PARAM("bcd-pos", &bcdpos, "String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of Spatial barcode sequences (in Read 1)")        
    LONG_STRING_PARAM("umi-pos", &umipos, "String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of UMI sequences (in Read 2)")
    //LONG_STRING_PARAM("tag-pos", &tagpos, "String in format of [beg1]-[end1],[beg2]-[end2],... to represent (1-based) positions of tag sequences (in Read 2)")

    LONG_PARAM_GROUP("Batch options", NULL)
    LONG_INT_PARAM("batch", &batch_size, "Size of single batch to store stored files")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outprefix,"Output prefix to store the tags")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( fq1f.empty() || fq2f.empty() || tagf.empty() || outprefix.empty() ) {
    error("Missing required options --fq1, --fq2, --tag, or --out");
  }

  // parse the positions of spatial barcodes, UMIs, and tags
  std::vector<uint64_t> bcd_begs, bcd_ends, umi_begs, umi_ends, tag_begs, tag_ends;
  int32_t bcd_len = 0, umi_len = 0, tag_len = 0;
  if (!bcdpos.empty()) {
    if ( !str2intervals(bcd_begs, bcd_ends, bcdpos.c_str()) ) {
      error("Cannot parse --bcd-pos argument %s", bcdpos.c_str());
    }
    for(int32_t i=0; i < (int32_t)bcd_begs.size(); ++i) {
      bcd_len += (bcd_ends[i] - bcd_begs[i] + 1);
    }
  }
  if (!umipos.empty()) {
    if ( !str2intervals(umi_begs, umi_ends, umipos.c_str()) ) {
      error("Cannot parse --umi-pos argument %s", umipos.c_str());
    }
    for(int32_t i=0; i < (int32_t)umi_begs.size(); ++i) {
      umi_len += (umi_ends[i] - umi_begs[i] + 1);
    }    
  }
  notice("Lengths of barcode = %d, UMI = %d, tag = %d", bcd_len, umi_len, tag_len);

  // read the tag reference file
  tsv_reader tr(tagf.c_str());
  std::vector<std::string> tag_ids;
  std::vector<std::string> tag_names;
  std::vector<std::string> tag_seqs;
  std::set<std::string> tag_id_set;
  char tag_seq[255];
  int32_t tag_ipos[255];

  for(int32_t i=0, k=0; i < (int32_t)tag_begs.size(); ++i) {
    for(int32_t j=tag_begs[i]; j <= tag_ends[i]; ++j) {
      tag_ipos[k++] = j - 1;
    }
  }
  while( tr.read_line() ) {
    if ( tr.nfields != 3 )
      error("Expects 3 columns consisting of [TAG_ID] [TAG_NAME] [TAG_SEQ]");
    const char* tag_id   = tr.str_field_at(0);
    const char* tag_name = tr.str_field_at(1);
    const char* seq      = tr.str_field_at(2);

    tag_ids.push_back(tag_id);
    tag_names.push_back(tag_name);
    tag_seqs.push_back(seq);

    if ( !tag_id_set.insert(tag_id).second )
      error("Duplicated Tag ID %s is observed", tag_id);
  }
  tr.close();

  int32_t n_tags = (int32_t)tag_ids.size();
  notice("Searching for the following tags:");
  for(int32_t i=0; i < n_tags; ++i) {
    notice("[%s]\t[%s]\t[%s]", tag_ids[i].c_str(), tag_names[i].c_str(), tag_seqs[i].c_str());
  }

  // load the smatch file if exists
  std::set<uint64_t> bcd_pass; // set of barcodes to be used as pass list
  int32_t bcd_nt5_len = bcd_len > MAX_NT5_LEN ? MAX_NT5_LEN : bcd_len;
  if ( !smatchf.empty() ) {
    notice("Loading smatch file %s", smatchf.c_str());
    tsv_reader tr(smatchf.c_str());
    while( tr.read_line() ) {
      uint64_t nt5 = seq2nt5(tr.str_field_at(0), bcd_nt5_len);
      bcd_pass.insert(nt5);
    }
    notice("Finished loading %zu barcodes from smatch file %s", bcd_pass.size(), smatchf.c_str());    
  }  

  // NOTE: here we could consider expanding seq2idx to allow mismatches

  notice("Finished loading %zu tags with length %d", tag_ids.size(), tag_len);

  // Read the pair of FASTQ files
  htsFile* hp1 = hts_open(fq1f.c_str(), "r");
  htsFile* hp2 = hts_open(fq2f.c_str(), "r");
  
  notice("Reading FASTQ files %s and %s", fq1f.c_str(), fq2f.c_str());
  
  int32_t lstr1, lstr2, lseq1, lseq2, ldummy1, ldummy2, lqual1, lqual2;  
  kstring_t str1; str1.l = str1.m = 0; str1.s = NULL;
  lstr1 = hts_getline(hp1, KS_SEP_LINE, &str1);
  kstring_t str2; str2.l = str2.m = 0; str2.s = NULL;
  lstr2 = hts_getline(hp2, KS_SEP_LINE, &str2);
  uint64_t bcd_nt5, nrecs = 0;
  char umi_seq[255], bcd_seq[255];
  int32_t umi_ipos[255], bcd_ipos[255];

  // determine umi_ipos and bcd_ipos

  for(int32_t i=0, k=0; i < (int32_t)umi_begs.size(); ++i) {
    for(int32_t j=umi_begs[i]; j <= umi_ends[i]; ++j) {
      umi_ipos[k++] = j - 1;
    }
  }

  for(int32_t i=0, k=0; i < (int32_t)bcd_begs.size(); ++i) {
    for(int32_t j=bcd_begs[i]; j <= bcd_ends[i]; ++j) {
      bcd_ipos[k++] = j - 1;
    }
  }

  // read the FASTQ1 and FASTQ2 and write the output data
  uint64_t cnts[4] = {0, 0, 0, 0};
  bool tag_match, bcd_match;
  htsFile* wbatch = NULL;
  uint64_t ibatch = 0, n_written = 0;
  while( lstr1 > 0 ) {
    if ( nrecs % 1000000 == 0 ) 
      notice("Processing %llu records from the FASTQ files, tb: %llu (%.5lf), tB: %llu (%.5lf), Tb: %llu (%.5lf), TB: %llu (%.5lf)", nrecs, cnts[0], (double)cnts[0]/(double)nrecs, cnts[1], (double)cnts[1]/(double)nrecs, cnts[2], (double)cnts[2]/(double)nrecs, cnts[3], (double)cnts[3]/(double)nrecs);

    // read the sequence reads for line 4N+1
    lseq1 = hts_getline(hp1, KS_SEP_LINE, &str1);
    lseq2 = hts_getline(hp2, KS_SEP_LINE, &str2);
    if ( lseq1 < bcd_len )
      error("Cannot parse sequence in FASTQ file %s at record=%llu. Read length is too short (%d)", fq1f.c_str(), nrecs, lseq1);
    if ( lseq2 < tag_len )
      error("Cannot parse sequence in FASTQ file %s at record=%llu. Read length is too short (%d)", fq2f.c_str(), nrecs, lseq2);

    for(int32_t i=0; i < bcd_len; ++i) bcd_seq[i] = str1.s[bcd_ipos[i]];
    if ( umi_len > 0 ) {
      for(int32_t i=0; i < umi_len; ++i) umi_seq[i] = str2.s[umi_ipos[i]];      
    }
    bcd_seq[bcd_len] = umi_seq[umi_len] = '\0';
    bcd_nt5 = seq2nt5(bcd_seq, bcd_nt5_len);

    // iterate all tag sequences and see if there is a match
    int32_t itag = -1;
    int32_t ipos = -1;
    for(int32_t i=0; i < n_tags; ++i) {
      const char* ret = strstr(str2.s, tag_seqs[i].c_str());
      if ( ret != NULL ) {
        if ( itag >= 0 ) {
          notice("Ignoring sequence matching to multiple tags : %s", str2.s);
          itag = n_tags;
        }
        else {
          itag = i;
          ipos = ret - str2.s;
        }
      }
    }

    tag_match = ( itag >= 0 ) && ( itag < n_tags );
    bcd_match = smatchf.empty() || ( bcd_pass.find(bcd_nt5) != bcd_pass.end() );

    if ( tag_match && bcd_match ) { // if both matches, write to the output file
      // create a new batch output
      if ( n_written % batch_size == 0 ) {
        notice("Writing batch %llu", ibatch);
        char filename[65535];
        if ( ibatch > 0 )
          hts_close(wbatch);
        sprintf(filename, "%s.batch.%llu.unsorted.tsv", outprefix.c_str(), static_cast<unsigned long long>(ibatch));
        wbatch = hts_open(filename, "w");
        ++ibatch;
      }
      // write [BARCODE] [TAG_ID] [UMI] [POS]for matching entries
      if ( umi_len == 0 ) {
        umi_seq[0] = '.';
        umi_seq[1] = '\0';
      }
      hprintf(wbatch, "%s\t%06d\t%s\t%d\n", bcd_seq, itag, umi_seq, ipos);
      ++n_written;
    }
    
    ++cnts[tag_match * 2 + bcd_match];

    ++nrecs;
    
    ldummy1 = hts_getline(hp1, KS_SEP_LINE, &str1);
    ldummy2 = hts_getline(hp2, KS_SEP_LINE, &str2);
    
    lqual1  = hts_getline(hp1, KS_SEP_LINE, &str1);
    lqual2  = hts_getline(hp2, KS_SEP_LINE, &str2);    

    lstr1 = hts_getline(hp1, KS_SEP_LINE, &str1);
    lstr2 = hts_getline(hp2, KS_SEP_LINE, &str2);        
  }
  hts_close(hp1);
  hts_close(hp2);
  hts_close(wbatch);

  notice("Finished processing %llu records from the FASTQ files, tb: %llu (%.5lf), tB: %llu (%.5lf), Tb: %llu (%.5lf), TB: %llu (%.5lf)", nrecs, cnts[0], (double)cnts[0]/(double)nrecs, cnts[1], (double)cnts[1]/(double)nrecs, cnts[2], (double)cnts[2]/(double)nrecs, cnts[3], (double)cnts[3]/(double)nrecs);

  htsFile* wtsv = hts_open((outprefix + ".unsorted.manifest.tsv").c_str(), "w");
  for(uint64_t i=0; i < ibatch; ++i) {
    hprintf(wtsv, "%s.batch.%llu.unsorted.tsv\n", outprefix.c_str(), i);
  }
  hts_close(wtsv);  

  notice("Analysis finished");
  
  return 0;
}
