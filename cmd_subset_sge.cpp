#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "file_utils.h"
#include "sge.h"
#include "polygon.h"
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
  std::string bcdwhtf;                        // name of barcode whitelist file
  std::string jsonf;                          // name of geojson file containing polygon (in um scale)
  int32_t xmin = 0;                           // minimum x coordinate
  int32_t xmax = INT_MAX;                     // maximum x coordinate
  int32_t ymin = 0;                           // minimum y coordinate
  int32_t ymax = INT_MAX;                     // maximum y coordinate
  double json_x_offset = 0.0;                 // x-offset to add to the geojson boundary
  double json_y_offset = 0.0;                 // y-offset to add to the geojson boundary
  double px_per_um = 26.67;                   // pixels per um
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
  LONG_STRING_PARAM("json", &jsonf, "Geojson file containing multiple polygons")
  LONG_DOUBLE_PARAM("json-x-offset", &json_x_offset, "X-offset to add to the geojson boundary")
  LONG_DOUBLE_PARAM("json-y-offset", &json_y_offset, "Y-offset to add to the geojson boundary")
  LONG_STRING_PARAM("whitelist", &bcdwhtf, "Barcode whitelist file path")
  LONG_DOUBLE_PARAM("px-per-um", &px_per_um, "Pixels/um scale (default: 26.67)")

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

  // if barcode whitelist exists, load it
  std::set<std::string> bcdwht;
  if (!bcdwhtf.empty())
  {
    notice("Loading barcode whitelist from %s", bcdwhtf.c_str());
    tsv_reader tr_wht(bcdwhtf.c_str());
    while (tr_wht.read_line())
    {
      bcdwht.insert(tr_wht.str_field_at(0));
    }
    notice("Finished loading %zu barcodes from the whitelist", bcdwht.size());
  }

  // if json file exists, load it
  std::vector<Polygon> polygons;
  if (!jsonf.empty())
  {
    notice("Loading polygons from %s", jsonf.c_str());
    int32_t npolygons = load_polygons_from_geojson(jsonf.c_str(), polygons);
    for(int32_t i=0; i < npolygons; ++i) {    
      polygons[i].add_offset(json_x_offset, json_y_offset);
    }
    notice("Finished loading %zu polygons", polygons.size());
  }

  sge_stream_reader ssr((sgedir + "/" + bcdf).c_str(), (sgedir + "/" + ftrf).c_str(), (sgedir + "/" + mtxf).c_str());

  // create the output directory first
  if (makePath(outdir))
  {
    notice("Successfully created the directory %s", outdir.c_str());
  }
  else
  {
    notice("The directory %s already exists", outdir.c_str());
  }
  sge_stream_writer ssw((outdir + "/" + bcdf).c_str(), (outdir + "/" + ftrf).c_str(), (outdir + "/" + mtxf).c_str());

  // read matrix
  bool bcd_pass = false;
  uint64_t npass = 0, nfail = 0;
  uint64_t out_xmin = UINT64_MAX, out_xmax = 0, out_ymin = UINT64_MAX, out_ymax = 0;
  while (ssr.read_mtx())
  {
    if (ssr.is_bcd_new)
    {
      bcd_pass = (ssr.cur_sbcd.px >= xmin && ssr.cur_sbcd.px <= xmax && ssr.cur_sbcd.py >= ymin && ssr.cur_sbcd.py <= ymax);
      if (bcd_pass && !bcdwht.empty())
      {
        // if barcode whitelist is not empty, check if the barcode is in the whitelist
        bcd_pass = (bcdwht.find(ssr.cur_sbcd.strid) != bcdwht.end());
      }
      if (bcd_pass && !jsonf.empty())
      {
        bcd_pass = polygons_contain_point(polygons, ssr.cur_sbcd.py/px_per_um, ssr.cur_sbcd.px/px_per_um);
      }
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
      ++npass;
      ssw.add_mtx(ssr.cur_iftr, ssr.cur_cnts);
    }
    else
    {
      ++nfail;
    }
    if ( ( npass + nfail ) % 10000000 == 0 )
    {
      notice("Processed %zu barcodes, %zu (%.5f) passed, %zu (%.5f) failed", npass + nfail, npass, (double)npass/(double)(npass + nfail), nfail, (double)nfail/(double)(npass + nfail));
    }
  }
  // ssw.flush_cur_sbcd(); // write out the last barcode

  // load features
  int32_t nftrs = ssr.load_features();
  ssw.ftr_cnts.resize(nftrs);
  for (int32_t i = 0; i < nftrs; ++i)
  {
    ssw.write_ftr(ssr.ftrs[i]->id.c_str(), ssr.ftrs[i]->name.c_str(), (uint32_t)(i + 1), ssw.ftr_cnts[i]);
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
