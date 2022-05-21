#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "file_utils.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

class tile_writer {
public:
  std::map<int32_t, std::map<int32_t, htsFile*> > tile_fhs;
  std::map<int32_t, htsFile*> lane_fhs;  
  std::map<int32_t, std::string> lane_filenames;   
  std::map<int32_t, std::map<int32_t, std::string> > tile_filenames;    
  htsFile* all_fh;
  uint64_t all_cnt;
  std::string all_filename;
  std::string rootdir;
  std::string filename;
  bool gzip;
  
  tile_writer(const char* _rootdir, const char* _filename, bool _gzip) : gzip(_gzip) {
    rootdir.assign(_rootdir);
    if ( rootdir[rootdir.size()-1] != '/' )
      rootdir += "/";
    filename.assign(_filename);
    all_filename = rootdir + filename;
    all_fh = hts_open(all_filename.c_str(), gzip ? "wz" : "w");
  }

  void close() {
    if ( all_fh != NULL ) {
      hts_close(all_fh);
      all_fh = NULL;
    }
    for(std::map<int32_t, htsFile*>::iterator it = lane_fhs.begin(); it != lane_fhs.end(); ++it) {
      hts_close(it->second);
    }
    lane_fhs.clear();
    
    for(std::map<int32_t, std::map<int32_t, htsFile*> >::iterator it = tile_fhs.begin(); it != tile_fhs.end(); ++it) {
      for(std::map<int32_t, htsFile*>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        hts_close(it2->second);
      }
    }
    tile_fhs.clear();
  }

  /*
  bool write_each(int32_t lane, int32_t tile, const std::string& s) {
    std::map<int32_t, htsFile*>::iterator it = tile_fhs[lane].find(tile);
    bool ret = false;
    htsFile* wf = NULL;
    if ( it == tile_fhs[lane].end() ) {
      // check if the directory exists
      std::string name;
      catprintf(name, "%s%d/%d", rootdir.c_str(), lane, tile);      
      ret = makePath(name); // make directory if it does not exist
      catprintf(name, "/%s", filename.c_str());
      wf = hts_open(name.c_str(), gzip ? "wz" : "w");      
      tile_fhs[lane][tile] = wf;
      tile_filenames[lane][tile] = name;
    } else {
      wf = it->second;
    }
    //++tile_cnts[lane][tile];
    //++all_cnt;
    // write the string twice
    hprint_str(wf, s);
    hprint_str(all_fh, s);
    return ret;
  }
  */

  void write_all(const std::string& s) {
    //++all_cnt;
    // write the string twice
    hprint_str(all_fh, s);
  }

  bool write_lane(int32_t lane, const std::string& s) {
    std::map<int32_t, htsFile*>::iterator it = lane_fhs.find(lane);
    bool ret = false;
    htsFile* wf = NULL;
    if ( it == lane_fhs.end() ) {
      std::string name;
      catprintf(name, "%s%d", rootdir.c_str(), lane);
      ret = makePath(name); // make directory if it does not exist
      catprintf(name, "/%s", filename.c_str());
      wf = hts_open(name.c_str(), gzip ? "wz" : "w");      
      lane_fhs[lane] = wf;
      lane_filenames[lane] = name;      
    } else {
      wf = it->second;
    }
    hprint_str(wf, s);
    return ret;    
  }

  bool write_tile(int32_t lane, int32_t tile, const std::string& s) {
    std::map<int32_t, htsFile*>::iterator it = tile_fhs[lane].find(tile);
    bool ret = false;
    htsFile* wf = NULL;
    if ( it == tile_fhs[lane].end() ) {
      std::string name;
      catprintf(name, "%s%d/%d", rootdir.c_str(), lane, tile);
      ret = makePath(name); // make directory if it does not exist
      catprintf(name, "/%s", filename.c_str());
      wf = hts_open(name.c_str(), gzip ? "wz" : "w");      
      tile_fhs[lane][tile] = wf;
      tile_filenames[lane][tile] = name;      
    } else {
      wf = it->second;
    }
    hprint_str(wf, s);
    return ret;
  }  
};

class tile_counter {
public:
  std::map<int32_t, std::map<int32_t, std::vector<uint64_t> > > tile_cnts;
  std::map<int32_t, std::vector<uint64_t> > lane_cnts;
  std::vector<uint64_t> all_cnts;
  int32_t ncols;
  
  tile_counter(int _ncols) : ncols(_ncols) {
    all_cnts.resize(ncols, 0);
  }

  inline uint64_t get_tile_count(int32_t lane, int32_t tile, int32_t icol) {
    std::vector<uint64_t>& cnts = tile_cnts[lane][tile];
    if ( cnts.empty() ) {
      cnts.resize(ncols,0);
    }
    return cnts[icol];
  }

  inline uint64_t get_lane_count(int32_t lane, int32_t icol) {
    std::vector<uint64_t>& cnts = lane_cnts[lane];
    if ( cnts.empty() ) {
      cnts.resize(ncols,0);
    }
    return cnts[icol];    
  }

  bool add_count(int32_t lane, int32_t tile, int32_t icol, int32_t cnt) {
    all_cnts[icol] += cnt;

    std::vector<uint64_t>& lcnts = lane_cnts[lane];
    if ( lcnts.empty() ) {
      lcnts.resize(ncols,0);
      lcnts[icol] = cnt;
    }
    else {
      lcnts[icol] += cnt;
    }    
    
    std::vector<uint64_t>& tcnts = tile_cnts[lane][tile];
    if ( tcnts.empty() ) {
      tcnts.resize(ncols,0);
      tcnts[icol] = cnt;
      return true;
    }
    else {
      tcnts[icol] += cnt;
      return false;
    }
  }

  bool add_counts(int32_t lane, int32_t tile, std::vector<int32_t>& v) {
    for(int32_t i=0; i < ncols; ++i)
      all_cnts[i] += v[i];

    std::vector<uint64_t>& lcnts = lane_cnts[lane];
    if ( lcnts.empty() ) {
      lcnts.resize(ncols,0);
      for(int32_t i=0; i < ncols; ++i)
        lcnts[i] = v[i];
    }
    else {
      for(int32_t i=0; i < ncols; ++i)
        lcnts[i] += v[i];      
    }    
            
    std::vector<uint64_t>& tcnts = tile_cnts[lane][tile];
    if ( tcnts.empty() ) {
      tcnts.resize(ncols,0);
      for(int32_t i=0; i < ncols; ++i)
        tcnts[i] = v[i];
      return true;
    }
    else {
      for(int32_t i=0; i < ncols; ++i)
        tcnts[i] += v[i];      
      return false;
    }
  }   
};

class sbcd_sync_reader {
public:
  std::vector<tsv_reader*> sbcd_trs;
  tsv_reader bcd_tr;
  int32_t idx_match;
  int32_t lbcd;
  const char* bcd;

  sbcd_sync_reader(const char* bcdfile) {
    bcd_tr.open(bcdfile);
    bcd_tr.read_line();
  }

  ~sbcd_sync_reader() { close(); }

  void close() {
    for(int32_t i=0; i < (int32_t)sbcd_trs.size(); ++i)
      sbcd_trs[i]->close();
    bcd_tr.close();
    sbcd_trs.clear();    
  }

  void add_sbcd_file(const char* sbcdfile) {
    tsv_reader *p = new tsv_reader(sbcdfile); // sorted barcode file
    p->read_line(); // peek each line
    sbcd_trs.push_back(p);
  }

  tsv_reader* move_to_ibcd(uint64_t ibcd) {
    // find the right barcode first
    while( bcd_tr.nfields > 0 && bcd_tr.nlines < ibcd ) {
      bcd_tr.read_line();
    }
    bcd = bcd_tr.str_field_at(0);
    lbcd = strlen(bcd);

    // scan the matching barcodes from the sbcd files
    int32_t nsbcds = (int32_t) sbcd_trs.size();
    int32_t cmp;
    idx_match = -1;
    for(int32_t i = 0; i < nsbcds; ++i) {
      if ( sbcd_trs[i]->nfields > 0 ) { // consider only when EOF is not yet reached
        while( ( cmp = strncmp(bcd, sbcd_trs[i]->str_field_at(0), lbcd) ) > 0 ) {
          if ( sbcd_trs[i]->read_line() == 0 ) {
            cmp = -1;
            break;
          }
        }
        if ( cmp == 0 ) {
          // match found
          idx_match = i;
          return sbcd_trs[idx_match];
          //return true;
        }
      }
    }
    return NULL;
    //return false; // match not found
  }
};

// once a single barcode finishes reading, flush the barcode out to all, lane, tile
void write_sbcd(sbcd_sync_reader& ssr, uint64_t cur_ibcd, std::vector<int32_t>& cur_bcd_cnts, tile_writer& bcd_tw, tile_counter& sbcds_counter, int32_t n_mtx) {
  // print the current barcode
  tsv_reader* p = ssr.move_to_ibcd(cur_ibcd);
  int32_t lane = p == NULL ? 0 : p->int_field_at(1);
  int32_t tile = p == NULL ? 0 : p->int_field_at(2);
  int32_t x    = p == NULL ? 0 : p->int_field_at(3);
  int32_t y    = p == NULL ? 0 : p->int_field_at(4);

  std::string outstr;

  // write all
  catprintf(outstr, "%s\t%llu\t%llu\t%d\t%d\t%d\t%d\t", ssr.bcd, sbcds_counter.all_cnts[0] + 1, cur_ibcd, lane, tile, x, y);  
  cat_join_int32(outstr, cur_bcd_cnts, ",");
  outstr += "\n";
  bcd_tw.write_all(outstr);

  outstr.clear();

  // write lane
  catprintf(outstr, "%s\t%llu\t%llu\t%d\t%d\t%d\t%d\t", ssr.bcd, sbcds_counter.get_lane_count(lane, 0) + 1, cur_ibcd, lane, tile, x, y);
  cat_join_int32(outstr, cur_bcd_cnts, ",");
  outstr += "\n";
  bcd_tw.write_lane(lane, outstr);

  outstr.clear();

  // write tile
  catprintf(outstr, "%s\t%llu\t%llu\t%d\t%d\t%d\t%d\t", ssr.bcd, sbcds_counter.get_tile_count(lane, tile, 0) + 1, cur_ibcd, lane, tile, x, y);  
  cat_join_int32(outstr, cur_bcd_cnts, ",");
  outstr += "\n";
  bcd_tw.write_tile(lane, tile, outstr);

  sbcds_counter.add_count(lane, tile, 0, 1); // add # barcodes
}

/////////////////////////////////////////////////////////////////////////
// dge2sdge : Convert multiple DGEs to a single SDGE
////////////////////////////////////////////////////////////////////////
int32_t cmdDGE2SDGE(int32_t argc, char** argv) {
  std::string sbcddir;
  std::vector<std::string> sbcdfs;
  std::string outdir;
  std::string bcdf; // barcodes.tsv
  std::string ftrf; // feature.tsv
  std::vector<std::string> mtxfs; // matrix.mtx files

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("sbcd", &sbcddir, "Spatial barcode dictionary generated from 'build-sbcds' command")
    LONG_MULTI_STRING_PARAM("match", &sbcdfs, "List of spatial barcode files that were used for whitelist generation")
    LONG_STRING_PARAM("bcd", &bcdf, "Shared barcode file path (e.g. barcodes.tsv.gz)")
    LONG_STRING_PARAM("ftr", &ftrf, "Shared feature file path (e.g. feature.tsv.gz)")
    LONG_MULTI_STRING_PARAM("mtx", &mtxfs, "Shared matrix file path (e.g. matrix.mtx.gz)")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outdir,"Output directory")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( bcdf.empty() || ftrf.empty() || mtxfs.empty() || outdir.empty() ) {
    error("Missing required options --bcd, --ftr, --mtx --out");
  }

  if ( !sbcddir.empty() ) {
    if ( !sbcdfs.empty() )
      error("--sbcd and --match cannot be used together");
        
    // make spatial barcode directory end with '/'
    if ( sbcddir[sbcddir.size()-1] != '/' )
      sbcddir += "/";

    // read the manifest file
    dataframe_t df((sbcddir + "manifest.tsv").c_str());
    // add a column to indicate the full path
    int32_t icol = df.add_empty_column("fullpath");
    int32_t jcol = df.get_colidx("filepath");
    
    for(int32_t i=0; i < df.nrows; ++i) {
      sbcdfs.push_back( sbcddir + df.get_str_elem(i, jcol).c_str() );
    }
  }
  else if ( sbcdfs.empty() )
    error("Either of --sbcd or --match parameters are required");

  // open spatial and regular barcodes
  sbcd_sync_reader ssr = (bcdf.c_str());
  int32_t nsbcds = (int32_t) sbcdfs.size();  
  for(int32_t i=0; i < nsbcds; ++i) {
    ssr.add_sbcd_file(sbcdfs[i].c_str());
  }

  int32_t n_mtx = (int32_t) mtxfs.size(); // number of input matrix files
  
  // open features.tsv files and set up feature counter
  tsv_reader ftr_tr(ftrf.c_str());
  std::vector<std::string> ftr_ids;
  std::vector<std::string> ftr_names;
  std::vector<tile_counter*> ftr_cnts;
  while( ftr_tr.read_line() > 0 ) {
    ftr_ids.push_back(ftr_tr.str_field_at(0));
    ftr_names.push_back(ftr_tr.str_field_at(1));
    ftr_cnts.push_back(new tile_counter(n_mtx));
  }

  // open all matrix.mtx files
  std::vector<tsv_reader*> mtx_trs;
  std::vector<bool> isValids;
  uint64_t nbcds = -1, nftrs = -1;
  for(int32_t i=0; i < (int32_t) mtxfs.size(); ++i) {   // scan each mtx file and read a single line
    tsv_reader* tr = new tsv_reader(mtxfs[i].c_str());
    mtx_trs.push_back(tr);
    while(tr->read_line() > 0) {
      if ( tr->str.s[0] == '%' ) {
          continue;
      }
      else {
        if ( i == 0 ) {
          nftrs = tr->uint64_field_at(0);
          nbcds = tr->uint64_field_at(1);
        }
        else {
          if ( nftrs != tr->uint64_field_at(0) )
            error("Number of features do not match %llu vs %llu", nftrs, tr->uint64_field_at(0));
          if ( nbcds != tr->uint64_field_at(1) )
            error("Number of features do not match %llu vs %llu", nbcds, tr->uint64_field_at(1));          
        }
        isValids.push_back( tr->read_line() > 0 );
        break;
      }
    }
  }

  // create the output directory first
  if ( makePath(outdir) ) {
    notice("Successfully created the directory %s", outdir.c_str());
  } else {
    notice("The directory %s already exists", outdir.c_str());    
  }
  
  // create a tile writer for matrix.mtx and barcode.tsv.gz
  tile_writer mtx_tw(outdir.c_str(), ".tmp.matrix.mtx", false);
  tile_writer bcd_tw(outdir.c_str(), "barcodes.tsv.gz", true);  

  // barcode cursors
  std::vector<int32_t> cur_bcd_cnts(n_mtx, 0); // counter for the current barcode

  //notice("foo1");

  // read over each matrix.mtx.gz files 
  bool hasValid = true, hasNZ = false;
  int32_t i, j, c;
  uint64_t ibcd, iftr, x, y;
  uint64_t cur_ibcd = 0;
  std::vector<int32_t> vals(n_mtx, 0);
  std::vector<bool> empty(n_mtx, false);
  tile_counter sbcds_counter(2); 
  std::string outstr;
  while( hasValid ) {
    hasNZ = false;
    ibcd = UINT64_MAX;
    iftr = UINT64_MAX;
    for(i=0; i < n_mtx; ++i) {
      if ( !isValids[i] ) {
        empty[i] = true;
        vals[i] = 0;
        continue;
      }
      x = mtx_trs[i]->uint64_field_at(0);
      y = mtx_trs[i]->uint64_field_at(1);
      c = mtx_trs[i]->int_field_at(2);
      if ( y < ibcd ) { // replace
        ibcd = y;
        iftr = x;
        for(j = 0; j < i; ++j) {
          vals[j] = 0;
          empty[j] = true;
        }
        vals[i] = c;
        empty[i] = false;
        hasNZ = c > 0;
      }
      else if ( y == ibcd ) {
        if ( x < iftr ) { // replace
          iftr = x;
          for(j = 0; j < i; ++j) {
            vals[j] = 0;
            empty[j] = true;
          }
          vals[i] = c;
          empty[i] = false;
          hasNZ = c > 0;
        }
        else if ( x == iftr ) { // add
          vals[i] = c;
          empty[i] = false;
          hasNZ |= (c > 0);
        }
        else {
          vals[i] = 0;
          empty[i] = true;
        }
      }
      else {
        vals[i] = 0;
        empty[i] = true;
      }
    }
    if ( hasNZ ) { // write output
      //notice("foo2");      
      if ( cur_ibcd < ibcd ) {
        if ( cur_ibcd > 0 ) { // flush out current barcode
          // print the current barcode
          write_sbcd(ssr, cur_ibcd, cur_bcd_cnts, bcd_tw, sbcds_counter, n_mtx);
          /*          tsv_reader* p = ssr.move_to_ibcd(cur_ibcd);
          int32_t lane = p == NULL ? 0 : p->int_field_at(1);
          int32_t tile = p == NULL ? 0 : p->int_field_at(2);
          int32_t x    = p == NULL ? 0 : p->int_field_at(3);
          int32_t y    = p == NULL ? 0 : p->int_field_at(4);
          int32_t sum  = cur_bcd_cnts[0];
          for(int32_t k=1; k < n_mtx; ++k) sum += cur_bcd_cnts[k];
          outstr.clear();
          catprintf(outstr, "%s\t%llu\t%llu\t%d\t%d\t%d\t%d\t%d\t", ssr.bcd, sbcds_written + 1, cur_ibcd, lane, tile, x, y, sum);
          cat_join_int32(outstr, cur_bcd_cnts, ",");
          outstr += "\n";
          bcd_tw.write_both(lane, tile, outstr);
          ++sbcds_written;
          */

          if ( sbcds_counter.all_cnts[0] % 500000 == 0 ) {
            notice("Processing %llu spatial barcodes, %llu/%llu = (%.5lf)", sbcds_counter.all_cnts[0], ibcd, nbcds, (double)ibcd/(double)nbcds);
          } 
        }
        cur_ibcd = ibcd;
        std::fill(cur_bcd_cnts.begin(), cur_bcd_cnts.end(), 0); // initialize the counts
      }
      outstr.clear();

      //notice("ibcd=%d, iftr=%d, cur_ibcd = %d", ibcd, iftr, cur_ibcd);
      tsv_reader* p = ssr.move_to_ibcd(ibcd);
      int32_t lane = p == NULL ? 0 : p->int_field_at(1);
      int32_t tile = p == NULL ? 0 : p->int_field_at(2);

      //notice("foo3");          

      // write matrix.mtx contents for merged file
      catprintf(outstr, "%llu %llu ", iftr, sbcds_counter.all_cnts[0] + 1);
      cat_join_int32(outstr, vals, " ");
      outstr += "\n";
      mtx_tw.write_all(outstr);
      outstr.clear();

      // write matrix.mtx contents for for specific lane
      catprintf(outstr, "%llu %llu ", iftr, sbcds_counter.get_lane_count(lane,0) + 1);
      cat_join_int32(outstr, vals, " ");
      outstr += "\n";
      mtx_tw.write_lane(lane, outstr);
      outstr.clear();      

      // write matrix.mtx contents for specific tile
      catprintf(outstr, "%llu %llu ", iftr, sbcds_counter.get_tile_count(lane,tile,0) + 1);
      cat_join_int32(outstr, vals, " ");
      outstr += "\n";
      mtx_tw.write_tile(lane, tile, outstr);
      outstr.clear();

      sbcds_counter.add_count(lane, tile, 1, 1); // add # lines

      //notice("foo5");    
      
      //mtx_tw.write_both(lane, tile, outstr);
      // update current barcode
      for(i=0; i < n_mtx; ++i) {
        cur_bcd_cnts[i] += vals[i];
      }
      // update current gene
      ftr_cnts[iftr-1]->add_counts(lane, tile, vals);
    }
    // advance the lines
    hasValid = false;
    for(i=0; i < n_mtx; ++i) {
      if ( !empty[i] ) {
        hasValid |= ( isValids[i] = mtx_trs[i]->read_line() > 0 );
      }
      else {
        hasValid |= isValids[i];
      }
    }
  }
  write_sbcd(ssr, cur_ibcd, cur_bcd_cnts, bcd_tw, sbcds_counter, n_mtx);
  
  ssr.close();
  bcd_tw.close();
  mtx_tw.close();
  
  notice("Finished wriing spatial barcodes and uncompressed .tmp.matrix.mtx files");

  notice("Started wriing compressed features.tsv.gz files");
  // write feature.tsv.gz file for each tile
  tile_writer ftr_tw( mtx_tw.rootdir.c_str(), "features.tsv.gz", true);
  for(int32_t i=0; i < (int32_t)ftr_ids.size(); ++i) {
    std::string buf1, buf2;
    tile_counter* pcnt = ftr_cnts[i];
    // buf1 contains basic gene info
    catprintf(buf1, "%s\t%s\t%d", ftr_ids[i].c_str(), ftr_names[i].c_str(), i+1);

    // buf2 contains actual counts
    catprintf(buf2, "%s\t", buf1.c_str());    
    cat_join_uint64(buf2, pcnt->all_cnts, ",");
    buf2 += "\n";
    ftr_tw.write_all(buf2); // write all file

    // write to each lane
    for(std::map<int32_t, std::vector<uint64_t> >::iterator it = pcnt->lane_cnts.begin(); it != pcnt->lane_cnts.end(); ++it) {
      buf2.clear();
      std::vector<uint64_t>& v = it->second;        
      catprintf(buf2, "%s\t", buf1.c_str());        
      cat_join_uint64(buf2, v, ",");
      buf2 += "\n";
      ftr_tw.write_lane(it->first, buf2);              
    }

    // write to each tile
    for(std::map<int32_t, std::map<int32_t, std::vector<uint64_t> > >::iterator it1 = pcnt->tile_cnts.begin(); it1 != pcnt->tile_cnts.end(); ++it1) {
      for(std::map<int32_t, std::vector<uint64_t> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
        buf2.clear();
        std::vector<uint64_t>& v = it2->second;        
        /* sum = 0;
        for(int32_t j=0; j < (int32_t)v.size(); ++j) {
          sum += v[j];
        }        
        catprintf(buf2, "%s\t%llu\t", buf1.c_str(), sum); */
        catprintf(buf2, "%s\t", buf1.c_str());        
        cat_join_uint64(buf2, v, ",");
        buf2 += "\n";
        ftr_tw.write_tile(it1->first, it2->first, buf2);        
      }
    }
  }
  ftr_tw.close();

  notice("Finished writing compressed features.tsv.gz files");

  // write header files
  tile_writer hdr_tw(outdir.c_str(), ".tmp.matrix.hdr", false);
  std::string hdr_str("%%MatrixMarket matrix coordinate integer general\n%\n");
  // write combined file
  outstr.clear();
  catprintf(outstr,"%s%d %llu %llu\n", hdr_str.c_str(), nftrs, sbcds_counter.all_cnts[0], sbcds_counter.all_cnts[1]);
  hdr_tw.write_all(outstr);
  // write each lane
  for(std::map<int32_t, std::vector<uint64_t> >::iterator it = sbcds_counter.lane_cnts.begin(); it != sbcds_counter.lane_cnts.end(); ++it) {
    outstr.clear();
    catprintf(outstr,"%s%d %llu %llu\n", hdr_str.c_str(), nftrs, sbcds_counter.get_lane_count(it->first, 0), sbcds_counter.get_lane_count(it->first, 1)); 
    hdr_tw.write_lane(it->first, outstr);      
  }  
  // write each tile
  for(std::map<int32_t, std::map<int32_t, std::vector<uint64_t> > >::iterator it1 = sbcds_counter.tile_cnts.begin(); it1 != sbcds_counter.tile_cnts.end(); ++it1) {
    for(std::map<int32_t, std::vector<uint64_t> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
      outstr.clear();
      catprintf(outstr,"%s%d %llu %llu\n", hdr_str.c_str(), nftrs, sbcds_counter.get_tile_count(it1->first, it2->first, 0), sbcds_counter.get_tile_count(it1->first, it2->first, 1)); 
      hdr_tw.write_tile(it1->first, it2->first, outstr);      
    }
  }
  hdr_tw.close();

  notice("Generating the merged matrix.mtx.gz file");  
  std::string cmd;
  catprintf(cmd, "cat %s %s | gzip -c > %s%s", hdr_tw.all_filename.c_str(), mtx_tw.all_filename.c_str(), mtx_tw.rootdir.c_str(), "matrix.mtx.gz");
  int32_t ret = system(cmd.c_str());
  if ( (ret == -1) || (WEXITSTATUS(ret) == 127) ) {
    error("Error in running %s", cmd.c_str());
  }
  if ( remove(hdr_tw.all_filename.c_str()) != 0 )
    error("Cannot remove %s", hdr_tw.all_filename.c_str());
  if ( remove(mtx_tw.all_filename.c_str()) != 0 )
    error("Cannot remove %s", mtx_tw.all_filename.c_str());

  notice("Generating the merged matrix.mtx.gz files for individual lanes");
  for(std::map<int32_t, std::vector<uint64_t> >::iterator it = sbcds_counter.lane_cnts.begin(); it != sbcds_counter.lane_cnts.end(); ++it) {
    cmd.clear();
    catprintf(cmd, "cat %s %s | gzip -c > %s%d/%s", hdr_tw.lane_filenames[it->first].c_str(), mtx_tw.lane_filenames[it->first].c_str(), mtx_tw.rootdir.c_str(), it->first, "matrix.mtx.gz");
    //notice("Lane %d, tile %d - %s", it1->first, it2->first, cmd.c_str());      
    ret = system(cmd.c_str());
    if ( (ret == -1) || (WEXITSTATUS(ret) == 127) ) {      
      error("Error in running %s", cmd.c_str());
    }
    if ( remove(hdr_tw.lane_filenames[it->first].c_str()) != 0 )
      error("Cannot remove %s", hdr_tw.lane_filenames[it->first].c_str());
    if ( remove(mtx_tw.lane_filenames[it->first].c_str()) != 0 )
      error("Cannot remove %s", mtx_tw.lane_filenames[it->first].c_str());        
  }    
  
  notice("Generating the merged matrix.mtx.gz files for individual tiles");
  for(std::map<int32_t, std::map<int32_t, std::vector<uint64_t> > >::iterator it1 = sbcds_counter.tile_cnts.begin(); it1 != sbcds_counter.tile_cnts.end(); ++it1) {
    for(std::map<int32_t, std::vector<uint64_t> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
      cmd.clear();
      catprintf(cmd, "cat %s %s | gzip -c > %s%d/%d/%s", hdr_tw.tile_filenames[it1->first][it2->first].c_str(), mtx_tw.tile_filenames[it1->first][it2->first].c_str(), mtx_tw.rootdir.c_str(), it1->first, it2->first, "matrix.mtx.gz");
      //notice("Lane %d, tile %d - %s", it1->first, it2->first, cmd.c_str());      
      ret = system(cmd.c_str());
      if ( (ret == -1) || (WEXITSTATUS(ret) == 127) ) {      
        error("Error in running %s", cmd.c_str());
      }
      if ( remove(hdr_tw.tile_filenames[it1->first][it2->first].c_str()) != 0 )
        error("Cannot remove %s", hdr_tw.tile_filenames[it1->first][it2->first].c_str());
      if ( remove(mtx_tw.tile_filenames[it1->first][it2->first].c_str()) != 0 )
        error("Cannot remove %s", mtx_tw.tile_filenames[it1->first][it2->first].c_str());        
    }
  }  
  
  notice("Analysis finished");
  
  return 0;
}
