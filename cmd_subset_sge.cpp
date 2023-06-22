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

/////////////////////////////////////////////////////////////////////////
// subset-sge : Subset individual SDGE file based on rectangular coordinates
////////////////////////////////////////////////////////////////////////
int32_t cmdSubsetSGE(int32_t argc, char **argv)
{
  std::string sgedir;                         // directory containg SGE file
  std::string outdir;                         // output directory
  std::string bcdf("barcodes.tsv.gz");        // name of barcode file
  std::string ftrf("features.tsv.gz");        // name of feature file
  std::string mtxf("matrix.mtx.gz");          // name of matrix file
  std::string minmaxf("barcodes.minmax.tsv"); // name of matrix file
  int32_t xmin = 0;                           // minimum x coordinate
  int32_t xmax = INT_MAX;                     // maximum x coordinate
  int32_t ymin = 0;                           // minimum y coordinate
  int32_t ymax = INT_MAX;                     // maximum y coordinate
  // TODO: later, we can add a whitelist option to subset based on barcode whitelist
  // TODO: later, we can add a boundarie polygon option to subset based on polygon

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
  LONG_PARAM_GROUP("Input options", NULL)
  LONG_STRING_PARAM("sge", &sgedir, "Spatial gene expression directory")
  LONG_STRING_PARAM("bcd", &bcdf, "Shared barcode file path (e.g. barcodes.tsv.gz)")
  LONG_STRING_PARAM("ftr", &ftrf, "Shared feature file path (e.g. feature.tsv.gz)")
  LONG_STRING_PARAM("mtx", &mtxf, "Shared matrix file path (e.g. matrix.mtx.gz)")

  LONG_PARAM_GROUP("Filter options", NULL)
  LONG_INT_PARAM("xmin", &xmin, "Minimum x coordinate")
  LONG_INT_PARAM("xmax", &xmax, "Maximum x coordinate")
  LONG_INT_PARAM("ymin", &ymin, "Minimum y coordinate")
  LONG_INT_PARAM("ymax", &ymax, "Maximum y coordinate")

  LONG_PARAM_GROUP("Output Options", NULL)
  LONG_STRING_PARAM("out", &outdir, "Output directory")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if (bcdf.empty() || ftrf.empty() || mtxf.empty() || outdir.empty())
  {
    error("Missing required options --bcd, --ftr, --mtx --out");
  }

  sge_stream_reader ssr((sgedir + "/" + bcdf).c_str(), (sgedir + "/" + ftrf).c_str(), (sgedir + "/" + mtxf).c_str());

  // create the output directory first
  if ( makePath(outdir) ) {
    notice("Successfully created the directory %s", outdir.c_str());
  } else {
    notice("The directory %s already exists", outdir.c_str());    
  }
  sge_stream_writer ssw((outdir + "/" + bcdf).c_str(), (outdir + "/" + ftrf).c_str(), (outdir + "/" + mtxf).c_str());

  // read matrix
  bool bcd_pass = false;
  uint64_t out_xmin = UINT64_MAX, out_xmax = 0, out_ymin = UINT64_MAX, out_ymax = 0;
  while (ssr.read_mtx())
  {
    if (ssr.is_bcd_new)
    {
      bcd_pass = (ssr.cur_sbcd.px >= xmin && ssr.cur_sbcd.px <= xmax && ssr.cur_sbcd.py >= ymin && ssr.cur_sbcd.py <= ymax);
      if (bcd_pass)
      {
        ssw.add_sbcd(ssr.cur_sbcd.strid.c_str(), ssr.cur_sbcd.nid, ssr.cur_sbcd.lane,
                     ssr.cur_sbcd.tile, ssr.cur_sbcd.px - xmin, ssr.cur_sbcd.py - ymin);
        // calculate the new bounding box
        if (ssr.cur_sbcd.px - (uint64_t)xmin > out_xmax)
          out_xmax = ssr.cur_sbcd.px - (uint64_t)xmin;
        if (ssr.cur_sbcd.px - (uint64_t)xmin < out_xmin)
          out_xmin = ssr.cur_sbcd.px - (uint64_t)xmin;
        if (ssr.cur_sbcd.py - (uint64_t)ymin > out_ymax)
          out_ymax = ssr.cur_sbcd.py - (uint64_t)ymin;
        if (ssr.cur_sbcd.py - (uint64_t)ymin < out_ymin)
          out_ymin = ssr.cur_sbcd.py - (uint64_t)ymin;
      }
    }
    if (bcd_pass)
    {
      ssw.add_mtx(ssr.cur_iftr, ssr.cur_cnts);
    }
  }
  //ssw.flush_cur_sbcd(); // write out the last barcode

  // load features
  int32_t nftrs = ssr.load_features();
  ssw.ftr_cnts.resize(nftrs);
  for (int32_t i = 0; i < nftrs; ++i)
  {
    ssw.write_ftr(ssr.ftrs[i].id.c_str(), ssr.ftrs[i].name.c_str(), (uint32_t)(i+1), ssw.ftr_cnts[i]);
  }

  ssw.close();

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
