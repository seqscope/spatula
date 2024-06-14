#ifndef __TILES_H__
#define __TILES_H__

#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "file_utils.h"
#include <ctime>
#include <set>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

class tile_writer
{
public:
  std::map<int32_t, std::map<int32_t, htsFile *>> tile_fhs;
  std::map<int32_t, htsFile *> lane_fhs;
  std::map<int32_t, std::string> lane_filenames;
  std::map<int32_t, std::map<int32_t, std::string>> tile_filenames;
  htsFile *all_fh;
  uint64_t all_cnt;
  std::string all_filename;
  std::string rootdir;
  std::string filename;
  bool gzip;

  tile_writer(const char *_rootdir, const char *_filename, bool _gzip);

  void close();

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

  inline void write_all(const std::string &s)
  {
    //++all_cnt;
    // write the string twice
    hprint_str(all_fh, s);
  }

  bool write_lane(int32_t lane, const std::string &s);

  bool write_tile(int32_t lane, int32_t tile, const std::string &s);
};

class tile_counter
{
public:
  std::map<int32_t, std::map<int32_t, std::vector<uint64_t>>> tile_cnts;
  std::map<int32_t, std::vector<uint64_t>> lane_cnts;
  std::vector<uint64_t> all_cnts;
  int32_t ncols;

  tile_counter(int _ncols) : ncols(_ncols)
  {
    all_cnts.resize(ncols, 0);
  }

  inline uint64_t get_tile_count(int32_t lane, int32_t tile, int32_t icol)
  {
    std::vector<uint64_t> &cnts = tile_cnts[lane][tile];
    if (cnts.empty())
    {
      cnts.resize(ncols, 0);
    }
    return cnts[icol];
  }

  inline uint64_t get_lane_count(int32_t lane, int32_t icol)
  {
    std::vector<uint64_t> &cnts = lane_cnts[lane];
    if (cnts.empty())
    {
      cnts.resize(ncols, 0);
    }
    return cnts[icol];
  }

  bool add_count(int32_t lane, int32_t tile, int32_t icol, int32_t cnt);

  bool add_counts(int32_t lane, int32_t tile, std::vector<int32_t> &v);
};

class sbcd_sync_reader
{
public:
  std::vector<tsv_reader *> sbcd_trs;
  tsv_reader bcd_tr;
  int32_t idx_match;
  int32_t lbcd;
  const char *bcd;

  sbcd_sync_reader(const char *bcdfile)
  {
    bcd_tr.open(bcdfile);
    bcd_tr.read_line();
  }

  ~sbcd_sync_reader() { close(); }

  void close()
  {
    for (int32_t i = 0; i < (int32_t)sbcd_trs.size(); ++i)
      sbcd_trs[i]->close();
    bcd_tr.close();
    sbcd_trs.clear();
  }

  void add_sbcd_file(const char *sbcdfile)
  {
    tsv_reader *p = new tsv_reader(sbcdfile); // sorted barcode file
    p->read_line();                           // peek each line
    sbcd_trs.push_back(p);
  }

  tsv_reader *move_to_ibcd(uint64_t ibcd);
};

struct _sbcd_rec_t {
  uint64_t nid;
  std::string strid;
  uint64_t lane;
  uint64_t tile;
  uint64_t px;
  uint64_t py;
  int32_t mismatches;

  _sbcd_rec_t(uint64_t _nid, const char* _strid, int32_t _lane, int32_t _tile, uint64_t _px, uint64_t _py, int32_t _mismatches) : 
    nid(_nid), strid(_strid), lane(_lane), tile(_tile), px(_px), py(_py), mismatches(_mismatches) {}

  void hprint_sbcd(htsFile* wh) {
    hprintf(wh, "%s\t%llu\t%llu\t%llu\t%llu\t%d\n", strid.c_str(), lane, tile, px, py, mismatches);
  }
};

typedef struct _sbcd_rec_t sbcd_rec_t;


void write_sbcd(sbcd_sync_reader &ssr, uint64_t cur_ibcd, std::vector<int32_t> &cur_bcd_cnts, tile_writer &bcd_tw, tile_counter &sbcds_counter, int32_t n_mtx);

std::pair<uint64_t,uint64_t> count_matches(std::vector<uint64_t>& bseqs, dataframe_t& df, std::vector<uint64_t>& dcounts, int32_t match_len, htsFile* wmatch);

uint64_t read_bcdf(tsv_reader* bcdf, int32_t match_len, uint64_t& cnt);
std::pair<uint64_t, uint64_t> count_matches_skip_dups(std::vector<uint64_t> &bseqs, dataframe_t &df, std::vector<uint64_t> &dcounts, int32_t match_len, htsFile *wmatch);

void open_tiles(dataframe_t& df, std::vector<std::string>& tiles, std::vector<tsv_reader*>& bcdfs);
void open_tiles(std::vector<std::string> &tile_paths, std::vector<tsv_reader*>& bcdfs);

#endif