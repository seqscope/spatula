#ifndef __SGE2_H__
#define __SGE2_H__

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

struct _sge2_sbcd_t
{
  std::string strid;  // string ID of barcode - usually the barcode sequence
  uint64_t nid;       // 1-based ID, automatically generated
  double px;          // spatial coordinates of X axis
  double py;          // spatial coordinate of Y axis
  std::vector<uint64_t> cnts; // counts of features

  _sge2_sbcd_t() : nid(0), px(0.0), py(0.0) {}

  void update_values(const char* _strid, double _px, double _py, uint64_t _nid = 0) {
    strid.assign(_strid);
    px = _px;
    py = _py;
    if ( _nid > 0 ) nid = _nid;
    else ++nid;

    // debug message - temporary
    if (rand() % 5000000 == 0)
      notice("strid = %s, nid = %llu, px = %llu, py = %llu", strid.c_str(), nid, px, py);
  }
};

typedef struct _sge2_sbcd_t sge2_sbcd_t;

struct _sge2_ftr_t
{
  std::string id;     // feature ID, must be unique
  std::string name;   // feature name, can be non-unique
  uint64_t nid;       // 1-based ID, automatically generated

  void update_values(const char* _id, const char* _name, uint64_t _nid = 0) {
    id.assign(_id);
    name.assign(_name);
    if ( _nid > 0 ) nid = _nid;
    else ++nid;
  }

  _sge2_ftr_t() : nid(0) {}

  _sge2_ftr_t(const char *_id, const char *_name, uint64_t _nid = 0) {
    update_values(_id, _name, _nid);
  }
};

typedef struct _sge2_ftr_t sge2_ftr_t;

// read SGE file in a streamed fashion
// matrix.mtx must be sorted by barcode first, and feature second
class sge2_stream_reader
{
public:
  tsv_reader bcd_tr;  // barcode file reader
  tsv_reader ftr_tr;  // feature file reader
  tsv_reader mtx_tr;  // matrix file reader
  std::map<std::string, std::pair<double, double> > bcd_pos_map;  // barcode position map (barcode -> (x, y))
  int32_t icol_bcd_strid;
  int32_t icol_bcd_px;
  int32_t icol_bcd_py;
  int32_t icol_ftr_id;
  int32_t icol_ftr_name;
  int32_t icol_mtx_cnts;

  // for barcodes.tsv.gz and features.tsv.gz, if the first line starts with #, it is considered as header
  // std::string colname_bcd_strid;
  // std::string colname_bcd_px;
  // std::string colname_bcd_py;
  // std::string colname_ftr_id;
  // std::string colname_ftr_name;

  int32_t nfields; // number of fields in the mtx file
  uint64_t nbcds;  // number of barcodes found in the mtx header
  uint64_t nftrs;  // number of features found in the mtx header
  uint64_t nlines; // number of lines found in the mtx hreader

  sge2_sbcd_t cur_sbcd; // current barcode cursor
  uint64_t cur_iftr;   // current index of feature in mtx cursor (1-based)
  uint64_t cur_line;   // current line number in mtx file
  bool is_bcd_new;
  std::vector<uint64_t> cur_cnts; // current counts in the mtx file
  std::vector<sge2_ftr_t*> ftrs;  // list of features

  void open(const char *bcdf, const char *ftrf, const char *mtxf);
  void close();

  sge2_stream_reader() {
    // header, if exists, must be the first line starting with #
    // colname_bcd_strid.assign("barcode");
    // colname_bcd_strid.assign("px");
    // colname_bcd_strid.assign("py");
    // colname_ftr_id.assign("id");
    // colname_ftr_name.assign("name");
    // default column indices, if header does not exists
    icol_bcd_strid = 0;
    icol_bcd_px = 1;
    icol_bcd_py = 2;
    icol_ftr_id = 0;
    icol_ftr_name = 1;
  }
  ~sge2_stream_reader() { close(); }

  bool read_mtx();
  int32_t load_features();
  uint64_t load_position_file(const char* filename, const char* colname_strid, const char* colname_px, const char* colname_py, char sep = '\t');
  uint64_t load_position_file(const char* filename, int32_t icol_strid, int32_t icol_px, int32_t icol_py, bool skip_header = false, char sep = '\t');
  uint64_t load_position_file_helper(tsv_reader& pos_tr, int32_t icol_strid, int32_t icol_px, int32_t icol_py, char sep = '\t'); // internal function
};

// read SGE file in a streamed fashion
// matrix.mtx must be provided in a sorted manner by barcode first, and feature second
class sge2_stream_writer
{
public:
  htsFile *wh_bcd;
  htsFile *wh_ftr;
  htsFile *wh_tmp;
  std::string fn_mtx;

  int32_t nfields;
  sge2_sbcd_t cur_sbcd;
  uint64_t nlines;
  int32_t bcd_precision;
  std::vector< std::vector<uint64_t> > ftr_cnts;

  void open(const char *bcdf, const char *ftrf, const char *mtxf);
  void close();

  sge2_stream_writer() { nfields = 0; bcd_precision = 3; nlines = 0; }
  sge2_stream_writer(const char *bcdf, const char *ftrf, const char *mtxf) { 
    nfields = 0; bcd_precision = 3; nlines = 0; 
    open(bcdf, ftrf, mtxf); 
  }
  ~sge2_stream_writer() { close(); }

  bool add_sbcd(const char *strid, double px, double py);
  bool add_mtx(uint64_t iftr, std::vector<uint64_t> &cnts);
  bool add_mtx(uint64_t iftr, std::vector<uint64_t> &cnts, std::vector<int32_t> &icols);
  bool flush_cur_sbcd();

  bool write_ftr(const char *id, const char *name, uint64_t nid, std::vector<uint64_t> &cnts);
  bool flush_mtx();
};

#endif