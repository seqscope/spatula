#ifndef __SGE_H__
#define __SGE_H__

#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/dataframe.h"
#include "seq_utils.h"
#include "file_utils.h"
#include "sge.h"
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

struct _sge_sbcd_t
{
  std::string strid;
  uint64_t nid;
  uint64_t gid;
  uint32_t lane;
  uint32_t tile;
  uint64_t px;
  uint64_t py;
  std::vector<uint64_t> cnts;

  _sge_sbcd_t() : nid(0), gid(0), lane(0), tile(0), px(0), py(0) {}

  void assign_fields(tsv_reader &tr)
  {
    strid.assign(tr.str_field_at(0));
    nid = tr.uint64_field_at(1);
    gid = tr.uint64_field_at(2);
    lane = (uint32_t)tr.int_field_at(3);
    tile = (uint32_t)tr.int_field_at(4);
    px = tr.uint64_field_at(5);
    py = tr.uint64_field_at(6);
    if (rand() % 5000000 == 0)
      notice("strid = %s, nid = %llu, px = %llu, py = %llu", strid.c_str(), nid, px, py);
  }
};

typedef struct _sge_sbcd_t sge_sbcd_t;

struct _sge_ftr_t
{
  std::string id;
  std::string name;
  uint64_t nid;
  std::vector<uint64_t> cnts;

  _sge_ftr_t() : nid(0) {}
  _sge_ftr_t(const char *_id, const char *_name, uint64_t _nid, const char *_cnts) : id(_id), name(_name), nid(_nid)
  {
    const char *pch = _cnts;
    while (pch != NULL)
    {
      uint64_t x = (uint64_t)strtoull(pch, NULL, 10);
      cnts.push_back(x);
      pch = strchr(pch, ',');
      if (pch != NULL)
        ++pch;
    }
  }
};

typedef struct _sge_ftr_t sge_ftr_t;

// read SGE file in a streamed fashion
// matrix.mtx must be sorted by barcode first, and feature second
class sge_stream_reader
{
public:
  tsv_reader bcd_tr;
  tsv_reader ftr_tr;
  tsv_reader mtx_tr;

  int32_t nfields; // number of fields in the mtx file
  uint64_t nbcds;  // number of barcodes in the mtx header
  uint64_t nftrs;  // number of features in the mtx header
  uint64_t nlines; // number of lines in the mtx hreader

  sge_sbcd_t cur_sbcd; // current barcode cursor
  uint64_t cur_iftr;   // current index of feature in mtx cursor
  uint64_t cur_line;   // current line number in mtx file
  bool is_bcd_new;
  std::vector<uint64_t> cur_cnts; // current counts in the mtx file
  std::vector<sge_ftr_t *> ftrs;  // list of features

  void open(const char *bcdf, const char *ftrf, const char *mtxf);
  void close();

  sge_stream_reader() {}
  sge_stream_reader(const char *bcdf, const char *ftrf, const char *mtxf) { open(bcdf, ftrf, mtxf); }
  ~sge_stream_reader() { close(); }

  bool read_mtx();
  int32_t load_features();
};

// read SGE file in a streamed fashion
// matrix.mtx must be provided in a sorted manner by barcode first, and feature second
class sge_stream_writer
{
public:
  htsFile *wh_bcd;
  htsFile *wh_ftr;
  htsFile *wh_tmp;
  std::string fn_mtx;

  int32_t nfields;
  sge_sbcd_t cur_sbcd;
  uint64_t nlines;
  std::vector<std::vector<uint64_t>> ftr_cnts;

  void open(const char *bcdf, const char *ftrf, const char *mtxf);
  void close();

  sge_stream_writer() {}
  sge_stream_writer(const char *bcdf, const char *ftrf, const char *mtxf) { open(bcdf, ftrf, mtxf); }
  ~sge_stream_writer() { close(); }

  bool add_sbcd(const char *strid, uint64_t old_nid, uint32_t lane, uint32_t tile, uint64_t x, uint64_t y);
  bool add_mtx(uint64_t iftr, std::vector<uint64_t> &cnts);
  bool flush_cur_sbcd();

  bool write_ftr(const char *id, const char *name, uint64_t nid, std::vector<uint64_t> &cnts);
  bool flush_mtx();
};

void write_sbcd(sbcd_sync_reader &ssr, uint64_t cur_ibcd, std::vector<int32_t> &cur_bcd_cnts, tile_writer &bcd_tw, tile_counter &sbcds_counter, int32_t n_mtx);

bool read_minmax(const char *fn, uint64_t &xmin, uint64_t &xmax, uint64_t &ymin, uint64_t &ymax);

std::pair<uint64_t,uint64_t> count_matches(std::vector<uint64_t>& bseqs, dataframe_t& df, std::vector<uint64_t>& dcounts, int32_t match_len, htsFile* wmatch);

uint64_t read_bcdf(tsv_reader* bcdf, int32_t match_len, uint64_t& cnt);
std::pair<uint64_t, uint64_t> count_matches_skip_dups(std::vector<uint64_t> &bseqs, dataframe_t &df, std::vector<uint64_t> &dcounts, int32_t match_len, htsFile *wmatch);

void open_tiles(dataframe_t& df, std::vector<std::string>& tiles, std::vector<tsv_reader*>& bcdfs);

#endif