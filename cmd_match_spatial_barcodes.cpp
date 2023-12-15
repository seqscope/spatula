#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "sge.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

/////////////////////////////////////////////////////////////////////////
// match-sbceds : Match spatial barcodes
////////////////////////////////////////////////////////////////////////
int32_t cmdMatchSpatialBarcodes(int32_t argc, char** argv) {
  std::string fastqf;
  std::string bcddir;
  std::string outprefix;
  int32_t batch_size = 300000000; // number of records to be searched for at once
  bool skip_duplicates = false;
  //int32_t hdmi_len = 32;
  int32_t match_len = 27;
  int32_t skip_sbcd = 0;
  int32_t nthreads = 1;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("fq", &fastqf, "FASTQ file read 1 containing 2nd-seq spatial barcode")
    LONG_STRING_PARAM("sbcd", &bcddir, "Spatial barcode dictionary generated from 'build-sbcds' command")
    LONG_INT_PARAM("batch", &batch_size, "Size of a single batch")
    LONG_INT_PARAM("skip-sbcd", &skip_sbcd, "Skip first bases of spatial barcode (Read 1)")
    LONG_INT_PARAM("match-len", &match_len, "Length of HDMI spatial barcodes to require perfect matches")
    LONG_PARAM("skip-duplicates", &skip_duplicates, "Skip duplicate barcodes that occurs multiple times")
    LONG_INT_PARAM("threads", &nthreads, "Number of threads")    
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outprefix,"Output prefix (index.tsv, matches.tsv.gz)")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( fastqf.empty() || bcddir.empty() || outprefix.empty() ) {
    error("Missing required options --fq, --sbcd, --out");
  }

  // make spatial barcode directory end with '/'
  if ( bcddir[bcddir.size()-1] != '/' ) 
    bcddir += "/";

  // read the manifest file
  dataframe_t df((bcddir + "manifest.tsv").c_str());
  if ( df.nrows == 0 )
    error("Empty dataframe %smanifest.tsv", bcddir.c_str());

  notice("Successfully read the manifest file, containing %d rows and %d columns", df.nrows, df.ncols);
  
  // add a column to indicate the full path
  int32_t icol = df.add_empty_column("fullpath");
  int32_t jcol = df.get_colidx("filepath");
  for(int32_t i=0; i < df.nrows; ++i) {
    df.set_str_elem((bcddir + df.get_str_elem(i, jcol)).c_str(), i, icol);
  }

  // write output files as plain files (per batch), and merge them later on
  std::vector<std::string> batch_filenames;
  int32_t cur_batch = 0;

  char buf[65536];
  snprintf(buf, 65536, "%s.match.batch.%d.tsv", outprefix.c_str(), cur_batch);
  htsFile* wmatch = hts_open(buf, "w");
  if ( wmatch == NULL )
    error("Cannot open %s for writing", buf);
  batch_filenames.push_back(buf);

  // Read the 2nd-seq FASTQ
  htsFile* hp = hts_open(fastqf.c_str(), "r");
  notice("Reading FASTQ file %s", fastqf.c_str());  
  std::vector<uint64_t> bseqs;      // store the 2nd-sequence into an array after conveting to 64bit integer
  //std::vector<uint64_t> ucounts, dcounts;
  std::vector<uint64_t> dcounts;  
  uint64_t nrecs = 0, ibatch = 0;
  //uint64_t nmissdups = 0;
  uint64_t nmiss = 0, ndups_approx = 0;
  int32_t lstr, lseq, ldummy, lqual;  
  kstring_t str; str.l = str.m = 0; str.s = NULL;
  lstr = hts_getline(hp, KS_SEP_LINE, &str);
  while( lstr > 0 ) {
    if ( nrecs % 10000000 == 0 ) 
      notice("Processing %d records from the FASTQ file %s", nrecs, fastqf.c_str());

    // read the sequence reads for line 4N+1
    lseq = hts_getline(hp, KS_SEP_LINE, &str);
    if ( lseq < skip_sbcd + match_len )
      error("Cannot parse Readname in FASTQ file %s at record=%llu. Read length is too short (%d)", fastqf.c_str(), nrecs, lseq);
    
    bseqs.push_back(seq2nt5(str.s + skip_sbcd, match_len));
    ldummy = hts_getline(hp, KS_SEP_LINE, &str);
    lqual  = hts_getline(hp, KS_SEP_LINE, &str);

    if ( ++nrecs % batch_size == 0 ) {
      notice("Processing batch %d of %d sequences", ++ibatch, batch_size);
      std::pair<uint64_t,uint64_t> nmissdups;
      if ( skip_duplicates ) {
        nmissdups = count_matches_skip_dups(bseqs, df, dcounts, match_len, wmatch);
      }
      else {
        nmissdups = count_matches(bseqs, df, dcounts, match_len, wmatch);
      }
      nmiss += nmissdups.first;
      ndups_approx += nmissdups.second;
      bseqs.clear();

      // open a new file
      hts_close(wmatch);
      cur_batch = nrecs / batch_size;
      snprintf(buf, 65536, "%s.match.batch.%d.tsv", outprefix.c_str(), cur_batch);
      wmatch = hts_open(buf, "w");
      batch_filenames.push_back(buf);      
      if ( wmatch == NULL )
        error("Cannot open %s.matches.tsv.gz for writing", outprefix.c_str());      
    }
    lstr = hts_getline(hp, KS_SEP_LINE, &str);    
  }
  if ( bseqs.size() > 0 ) {
    notice("Processing the last batch %d containing %zu sequences", ++ibatch, bseqs.size());
    std::pair<uint64_t,uint64_t> nmissdups;
    if ( skip_duplicates ) {
      nmissdups = count_matches_skip_dups(bseqs, df, dcounts, match_len, wmatch);
    }
    else {
      nmissdups = count_matches(bseqs, df, dcounts, match_len, wmatch);
    }
    nmiss += nmissdups.first;
    ndups_approx += nmissdups.second;    
    bseqs.clear();
  }
  hts_close(wmatch);

  // merge across batches and identify unique sequences
  std::vector<tsv_reader*> batch_trs;
  int32_t nbatches = (int32_t)batch_filenames.size();
  std::vector<int32_t> cmps(nbatches, 0);
  int32_t ndups = 0, nEOFs = 0;
  for(int32_t i=0; i < nbatches; ++i) {
    batch_trs.push_back(new tsv_reader(batch_filenames[i].c_str()));
    if ( batch_trs.back()->read_line() == 0 ) {
        cmps[i] = 999;
        ++nEOFs;
    }
  }

  htsFile* wmerged = hts_open((outprefix+".match.sorted.uniq.tsv.gz").c_str(), "wz");
  std::vector<int32_t> vals(5);
  std::map< int32_t, std::map<int32_t, uint64_t> > uniq_cnts;
  int32_t i, j, k;
  // find the new minimum
  strcpy(buf, "zzz"); // fill in dummy max

  for(i=0; i < nbatches; ++i) {
    if ( cmps[i] == 999 ) // if EOF was reached, skip
        continue;
    cmps[i] = strcmp(batch_trs[i]->str_field_at(0), buf);
    if ( cmps[i] < 0 ) { // new minimum found
      for(j=0; j < i; ++j)
        if ( cmps[j] == 0 )
          cmps[j] = 1;
      cmps[i] = 0;
      strcpy(buf, batch_trs[i]->str_field_at(0));
      for(k=0; k < 5; ++k) vals[k] = batch_trs[i]->int_field_at(k+1);      
    }
  }

  while( nEOFs < nbatches ) {
    for(i=0; i < nbatches; ++i) {
      while( cmps[i] == 0 ) {  // keep counting duplicates
        ++ndups;
        if ( batch_trs[i]->read_line() > 0 ) { // next line exists
          cmps[i] = strcmp(batch_trs[i]->str_field_at(0), buf);
        }
        else { // next line does not exists
          cmps[i] = 999; // EOF reached
          ++nEOFs;
        }
      }
    }
    // print the output
    hprintf(wmerged, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", buf, vals[0], vals[1], vals[2], vals[3], vals[4], ndups);
    ++uniq_cnts[vals[0]][vals[1]];
    ndups = 0;

    // determine the next minimum
    buf[0] = '\0';
    for(i=0; i < nbatches; ++i) {
      if ( cmps[i] != 999 ) { // not reached EOF yet
        cmps[i] = strcmp(batch_trs[i]->str_field_at(0), buf);
        if ( buf[0] == '\0' || cmps[i] < 0 ) { // new minimum found
          for(j=0; j < i; ++j) {
            if ( cmps[j] == 0 ) {
              cmps[j] = 1;
            }
          }
          cmps[i] = 0;
          strcpy(buf, batch_trs[i]->str_field_at(0));
          for(k=0; k < 5; ++k) vals[k] = batch_trs[i]->int_field_at(k+1);
        }
      }
    }  
  }
  hts_close(wmerged);  
  // make a map between ID and lane/tile
  std::map<std::string, std::pair<int32_t,int32_t> > id2lanetile;
  for(std::map<int32_t, std::map<int32_t,uint64_t> >::iterator it1=uniq_cnts.begin(); it1 != uniq_cnts.end(); ++it1) {
    for(std::map<int32_t,uint64_t>::iterator it2=it1->second.begin(); it2 != it1->second.end(); ++it2) {
      snprintf(buf,65536, "%d_%d", it1->first, it2->first);
      id2lanetile[buf] = std::make_pair(it1->first, it2->first);
    }
  }

  // read manifest files
  htsFile* wf = hts_open((outprefix + ".counts.tsv").c_str(), "w");
  //hprintf(wf, "id\tfilepath\tbarcodes\tmatches\tmatches_uniq_per_batch\n");
  hprintf(wf, "id\tfilepath\tbarcodes\tmatches\tunique\n");
  uint64_t nuniq = 0;
  for( int32_t i=0; i < df.nrows; ++i) {
    std::pair<int32_t,int32_t> pair = id2lanetile[df.get_str_elem(i, "id")];
    hprintf(wf, "%s", df.get_str_elem(i, "id").c_str());
    hprintf(wf, "\t%s", df.get_str_elem(i, "filepath").c_str());
    hprintf(wf, "\t%llu", df.get_uint64_elem(i, "barcodes"));
    hprintf(wf, "\t%llu\t%llu\n", dcounts[i], uniq_cnts[pair.first][pair.second]);
    nuniq += uniq_cnts[pair.first][pair.second];
  }
  hts_close(wf);

  wf = hts_open((outprefix + ".summary.tsv").c_str(), "w");
  hprintf(wf, "Type\tReads\tFraction\n");
  hprintf(wf, "Total\t%llu\t%.5lf\n", nrecs, 1.0);
  hprintf(wf, "Miss\t%llu\t%.5lf\n",  nmiss, (double)nmiss/(double)nrecs);
  hprintf(wf, "Match\t%llu\t%.5lf\n",  nrecs-nmiss, 1.0-(double)nmiss/(double)nrecs);
//  hprintf(wf, "Dup(Approx)\t%llu\t%.5lf\n",  ndups_approx, (double)ndups_approx/(double)nrecs);
  hprintf(wf, "Unique\t%llu\t%.5lf\n",  nuniq, (double)nuniq/(double)nrecs);
  hprintf(wf, "Dup(Exact)\t%llu\t%.5lf\n",  nrecs-nmiss-nuniq, (double)(nrecs-nmiss-nuniq)/(double)nrecs);  
  hts_close(wf);
  
  notice("Total = %llu, Missed = %llu (%.5f), Dups(approx) = %llu (%.5f), Unique = %llu (%.5lf)", nrecs, nmiss, (double)nmiss/(double)nrecs,  ndups_approx, (double)ndups_approx/(double)nrecs, nuniq, (double)nuniq/(double)nrecs);

  for(i=0; i < nbatches; ++i) {
    batch_trs[i]->close();
    if ( remove(batch_filenames[i].c_str()) == 0 ) {
      notice("Successfully removed %s", batch_filenames[i].c_str());        
    }
    else {
      error("ERROR in removing file %s", batch_filenames[i].c_str());              
    }
  }
  notice("Finished writing the merged file");
  notice("Analysis finished");
  
  return 0;
}
