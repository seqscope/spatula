#ifndef __SGE_H__
#define __SGE_H__

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

struct _sge_sbcd_t
{
  std::string strid;
  uint64_t nid; // 1-based ID
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
  uint64_t cur_iftr;   // current index of feature in mtx cursor (1-based)
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

bool read_minmax(const char *fn, uint64_t &xmin, uint64_t &xmax, uint64_t &ymin, uint64_t &ymax);
bool read_minmax_double(const char *fn, double &xmin, double &xmax, double &ymin, double &ymax);


#endif