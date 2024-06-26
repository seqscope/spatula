#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "file_utils.h"
#include "sge.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

///////////////////////////////////////////////////////////////////////////////////////////////
// combine-sge : Combine multiple SDGE tiles with global coordinates based on specified layouts
///////////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdCombineSGE(int32_t argc, char **argv)
{
  std::string layoutf;                        // Layout file, each containing [sgedir] and [row]/[col] or [x_offset]/[y_offset] as columns
  std::string outdir;                         // output directory
  std::string bcdf("barcodes.tsv.gz");        // name of barcode file
  std::string ftrf("features.tsv.gz");        // name of feature file
  std::string mtxf("matrix.mtx.gz");          // name of matrix file
  std::string minmaxf("barcodes.minmax.tsv"); // name of minmax TSV file indicating the boundary
  bool out_minmax_fixed = false;              // do not update output minmax coordinates based on the observed points
  bool mode_rowcol = false;                   // mode row_col is false by default

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
  LONG_PARAM_GROUP("Input options", NULL)
  LONG_STRING_PARAM("layout", &layoutf, "Layout file, each containing [sgedir] and [row]/[col] or [x_offset]/[y_offset] as columns")
  LONG_STRING_PARAM("bcd", &bcdf, "Shared barcode file path (e.g. barcodes.tsv.gz)")
  LONG_STRING_PARAM("ftr", &ftrf, "Shared feature file path (e.g. feature.tsv.gz)")
  LONG_STRING_PARAM("mtx", &mtxf, "Shared matrix file path (e.g. matrix.mtx.gz)")
  LONG_STRING_PARAM("minmax", &minmaxf, "Shared minmax.tsv file path (e.g. barcodes.minmax.tsv) - required in [row]/[col] mode")
  LONG_PARAM("out-minmax-fixed", &out_minmax_fixed, "Do not update output minmax coordinates based on the observed points")
  LONG_PARAM("mode-rowcol", &mode_rowcol, "Activate rowcol mode (default: false)")

  LONG_PARAM_GROUP("Output Options", NULL)
  LONG_STRING_PARAM("out", &outdir, "Output directory")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if (bcdf.empty() || ftrf.empty() || mtxf.empty() || outdir.empty() || layoutf.empty() ) 
  {
    error("Missing required options --bcd, --ftr, --mtx, --out, --layout");
  }  

  // determine the mode
  dataframe_t df(layoutf.c_str());
  //bool mode_rowcol = false;
  bool layout_has_minmax = false;

  if ( mode_rowcol ) {
    if ( !( df.has_column("row") && df.has_column("col") ) ) {
      error("The layout file must contain [row]/[col] when --mode-rowcol is set");
    }
  }
  else {
    if ( ! (df.has_column("x_offset") && df.has_column("y_offset") ) ) {
      error("The layout file must contain [x_offset]/[y_offset] when --mode-rowcol is not set");
    }
  }
  // if ( df.has_column("row") && df.has_column("col") ) {
  //   mode_rowcol = true;
  // }
  // else if ( df.has_column("x_offset") && df.has_column("y_offset") ) {
  //   mode_rowcol = false;
  // }
  // else {
  //   error("The layout file must contain either [row]/[col] or [x_offset]/[y_offset] as columns");
  // }

  if ( df.has_column("xmin") && df.has_column("xmax") && df.has_column("ymin") && df.has_column("ymax") ) {
    layout_has_minmax = true;
  }

  std::vector<std::string> sgedirs;
  std::vector<uint64_t> x_offsets, y_offsets, xmins, xmaxs, ymins, ymaxs;

  // store the bounding box information for each tile
  std::vector<std::string>& col_sgedir = df.get_column("sgedir");
  int32_t i_xmin, i_xmax, i_ymin, i_ymax;
  if ( layout_has_minmax ) {
    i_xmin = df.get_colidx("xmin");
    i_xmax = df.get_colidx("xmax");
    i_ymin = df.get_colidx("ymin");
    i_ymax = df.get_colidx("ymax");
  }

  for(int32_t i=0; i < df.nrows; ++i) {
    uint64_t xmin, xmax, ymin, ymax;
    if ( layout_has_minmax ) {    // if the layout file has xmin, xmax, ymin, ymax, use them as bounding box
      xmin = df.get_uint64_elem(i, i_xmin);
      xmax = df.get_uint64_elem(i, i_xmax);
      ymin = df.get_uint64_elem(i, i_ymin);
      ymax = df.get_uint64_elem(i, i_ymax);
    }
    else { // otherwise, read the minmax file from the sgedir
      read_minmax((col_sgedir[i] + "/" + minmaxf).c_str(), xmin, xmax, ymin, ymax);
    }
    sgedirs.push_back(col_sgedir[i]);
    xmins.push_back(xmin);
    xmaxs.push_back(xmax);
    ymins.push_back(ymin);
    ymaxs.push_back(ymax);
  }

  // store the offsets for each tile
  uint64_t out_max_row_x = 0, out_max_col_y = 0;
  uint64_t out_xmin = UINT64_MAX, out_xmax = 0, out_ymin = UINT64_MAX, out_ymax = 0;
  if ( mode_rowcol ) {
    // calculate the offsets for each row and column
    std::vector<std::string>& col_row = df.get_column("row");
    std::vector<std::string>& col_col = df.get_column("col");
    std::vector<int32_t> rows, cols;
    
    // read row, col information
    int32_t max_row = 0, max_col = 0;
    for(int32_t i=0; i < df.nrows; ++i) {
      int32_t row = atoi(col_row[i].c_str()) - 1;
      int32_t col = atoi(col_col[i].c_str()) - 1;
      if ( row < 0 || col < 0 ) {
        error("Row and column indices must be positive integers in %s", layoutf.c_str());
      }
      rows.push_back(row);
      cols.push_back(col);
      if ( row > max_row ) max_row = row;
      if ( col > max_col ) max_col = col;
    } 
  
    // determine maximum row height and column widths
    std::vector<uint64_t> max_row_heights(max_row+1, 0);
    std::vector<uint64_t> max_col_widths(max_col+1, 0);
    for(int32_t i=0; i < df.nrows; ++i) {
      int32_t row = rows[i];
      int32_t col = cols[i];
      uint64_t row_height = xmaxs[i] - xmins[i] + 1;
      uint64_t col_width = ymaxs[i] - ymins[i] + 1;
      if ( row_height > max_row_heights[row] ) max_row_heights[row] = row_height;
      if ( col_width > max_col_widths[col] ) max_col_widths[col] = col_width;
    }

    std::vector<uint64_t> row_cumsum(max_row+2, 0);
    std::vector<uint64_t> col_cumsum(max_col+2, 0);
    for(int32_t i=0; i < max_row; ++i) {
      row_cumsum[i+1] = row_cumsum[i] + max_row_heights[i];
    }
    for(int32_t i=0; i < max_col; ++i) {
      col_cumsum[i+1] = col_cumsum[i] + max_col_widths[i];
    }

    uint64_t min_x_offset = UINT64_MAX, min_y_offset = UINT64_MAX;

    // calculte x_offset and y_offset for each tile
    for(int32_t i=0; i < df.nrows; ++i) {
      int32_t row = rows[i];
      int32_t col = cols[i];
      uint64_t x_offset = row_cumsum[row];
      uint64_t y_offset = col_cumsum[col];
      x_offsets.push_back(x_offset);
      y_offsets.push_back(y_offset);

      if ( x_offsets[i] < min_x_offset ) min_x_offset = x_offsets[i];
      if ( y_offsets[i] < min_y_offset ) min_y_offset = y_offsets[i];

      notice("row = %d, col = %d, height = %llu, width = %llu, x_offset = %llu, y_offset = %llu", 
              row, col, xmaxs[i] - xmins[i] + 1, ymaxs[i] - ymins[i] + 1, x_offset, y_offset);
    }

    out_max_row_x = row_cumsum[row_cumsum.size()-1];
    out_max_col_y = col_cumsum[col_cumsum.size()-1];

    out_xmin = min_x_offset;
    out_xmax = out_max_row_x;
    out_ymin = min_y_offset;
    out_ymax = out_max_col_y;
  }
  else {
    // fill in the offsets
    std::vector<std::string>& col_sgedir = df.get_column("sgedir");
    std::vector<std::string>& col_x_offset = df.get_column("x_offset");
    std::vector<std::string>& col_y_offset = df.get_column("y_offset");

    uint64_t min_x_offset = UINT64_MAX, min_y_offset = UINT64_MAX;
    uint64_t max_x_value = 0, max_y_value = 0;

    for (int32_t i=0; i < df.nrows; ++i) {
      sgedirs.push_back(col_sgedir[i]);
      x_offsets.push_back(std::stoull(col_x_offset[i].c_str(), NULL, 10));
      y_offsets.push_back(std::stoull(col_y_offset[i].c_str(), NULL, 10));

      if ( x_offsets[i] < min_x_offset ) min_x_offset = x_offsets[i];
      if ( y_offsets[i] < min_y_offset ) min_y_offset = y_offsets[i];
      if ( xmaxs[i] - xmins[i] + 1 + x_offsets[i] > max_x_value ) max_x_value = xmaxs[i] - xmins[i] + 1 + x_offsets[i];
      if ( ymaxs[i] - ymins[i] + 1 + y_offsets[i] > max_y_value ) max_y_value = ymaxs[i] - ymins[i] + 1 + y_offsets[i];
    }
    if ( out_minmax_fixed ) {
      //error("--out-minmax-fixed must be used with input files with row/col");
      out_xmin = min_x_offset;
      out_xmax = max_x_value;
      out_ymin = min_y_offset;
      out_ymax = max_y_value;
    }
  }

  notice("Finished calculating offsets for %d tiles", df.nrows);

  // create the output directory first
  if (makePath(outdir))
  {
    notice("Successfully created the directory %s", outdir.c_str());
  }
  else
  {
    notice("The directory %s already exists", outdir.c_str());
  }

  // read each tile and write to the output
  sge_stream_writer ssw((outdir + "/" + bcdf).c_str(), (outdir + "/" + ftrf).c_str(), (outdir + "/" + mtxf).c_str());
  // if ( out_minmax_fixed ) { // use the maximum range possible as minmax values
  //   out_xmin = 0;
  //   out_xmax = out_max_row_x;
  //   out_ymin = 0;
  //   out_ymax = out_max_col_y;
  // }
  int32_t nftrs = 0;
  // TODO: make sure that feature.tsv.gz are compatible. If not, need to find a way to merge them.
  char buf[65535];
  for(int32_t i=0; i < df.nrows; ++i) {
    notice("Processing SGE directory %s ...", sgedirs[i].c_str());
    //sge_stream_reader ssr((sgedirs[i] + "/" + bcdf).c_str(), (sgedirs[i] + "/" + ftrf).c_str(), (sgedirs[i] + "/" + mtxf).c_str());
    sge_stream_reader* pssr = new sge_stream_reader((sgedirs[i] + "/" + bcdf).c_str(), (sgedirs[i] + "/" + ftrf).c_str(), (sgedirs[i] + "/" + mtxf).c_str());
    sge_stream_reader& ssr = *pssr;
    while( ssr.read_mtx() ) {
      // remove lane and tile info and create global coordinates
      if ( ssr.is_bcd_new ) {
        uint64_t gx = ssr.cur_sbcd.px + x_offsets[i];
        uint64_t gy = ssr.cur_sbcd.py + y_offsets[i];
        gx = gx < xmins[i] ? 0 : gx - xmins[i];
        gy = gy < ymins[i] ? 0 : gy - ymins[i];
        snprintf(buf, 65535, "%d_%s", i+1, ssr.cur_sbcd.strid.c_str()); // create unique spatial barcode id
        ssw.add_sbcd(buf, ssr.cur_sbcd.nid, 1, 1, gx, gy);
        if ( gx < out_xmin ) out_xmin = gx;
        if ( gx > out_xmax ) out_xmax = gx;
        if ( gy < out_ymin ) out_ymin = gy;
        if ( gy > out_ymax ) out_ymax = gy;
      }
      ssw.add_mtx(ssr.cur_iftr, ssr.cur_cnts);
    }
    if ( nftrs == 0 ) {
      nftrs = ssr.load_features();
      ssw.ftr_cnts.resize(nftrs);
    }
    else {
      int32_t nftrs2 = ssr.load_features();
      if ( nftrs != nftrs2 ) {
        error("The number of features in %s is different from the previous tiles", (sgedirs[i] + "/" + ftrf).c_str());
      }
    }
    if ( i == df.nrows - 1 ) {
      for (int32_t j = 0; j < nftrs; ++j)
      {
        ssw.write_ftr(ssr.ftrs[j]->id.c_str(), ssr.ftrs[j]->name.c_str(), (uint64_t)(j + 1), ssw.ftr_cnts[j]);
      }
    }
    delete pssr;  
    //ssr.close();
  }
  ssw.close();

  notice("Finished writing the combine SGE files");

  // write out the minmax file
  htsFile *fh_minmax = hts_open((outdir + "/" + minmaxf).c_str(), "w");
  if (fh_minmax == NULL)
  {
    error("Cannot open file %s for writing", (outdir + "/" + minmaxf).c_str());
  }
  hprintf(fh_minmax, "xmin\txmax\tymin\tymax\n");
  hprintf(fh_minmax, "%llu\t%llu\t%llu\t%llu\n", out_xmin, out_xmax, out_ymin, out_ymax);
  hts_close(fh_minmax);

  notice("Analysis finished");

  return 0;
}
