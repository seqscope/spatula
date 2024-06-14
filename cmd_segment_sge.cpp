/////////////////////////////////////////////////////////////////////////
// WARNING: This code is experimental and not finished yet
////////////////////////////////////////////////////////////////////////

#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge2.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include "libnpy/npy.hpp"

/////////////////////////////////////////////////////////////////////////
// convert-sge : Convert SGE format to a TSV format
////////////////////////////////////////////////////////////////////////
int32_t cmdSegmentSGE(int32_t argc, char **argv)
{
    // core input files
    std::string in_sgedir;
    std::string in_npy;
    std::string out_sgedir;

    std::string bcdf("barcodes.tsv.gz");
    std::string ftrf("features.tsv.gz");
    std::string mtxf("matrix.mtx.gz");

    int32_t in_icol_bcd_barcode = 1;
    int32_t in_icol_bcd_px = 6; // 1-based column index in the barcode file to use as the X coordinate
    int32_t in_icol_bcd_py = 7; // 1-based column index in the barcode file to use as the Y coordinate
    int32_t in_icol_ftr_id = 1;
    int32_t in_icol_ftr_name = 2;
    int32_t in_icol_mtx = 1; // 1-based column index in the matrix file to use as the count
    
    double units_per_px = 1.0;  // Output conversion factor 

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("in-sge", &in_sgedir, "Input SGE directory")
    LONG_STRING_PARAM("in-npy", &in_npy, "Input SGE npy files")
    LONG_STRING_PARAM("out-sge", &out_sgedir, "Prefix of output SGE directory")

    LONG_PARAM_GROUP("Input Options for SGE input", NULL)
    LONG_STRING_PARAM("bcd", &bcdf, "Barcode file name")
    LONG_STRING_PARAM("ftr", &ftrf, "Feature file name")
    LONG_STRING_PARAM("mtx", &mtxf, "Matrix file name")
    LONG_INT_PARAM("icol-mtx", &in_icol_mtx, "1-based column index in the matrix file to use as the count")
    LONG_INT_PARAM("icol-bcd-barcode", &in_icol_bcd_barcode, "1-based column index of barcode in the barcode file")
    LONG_INT_PARAM("icol-bcd-x", &in_icol_bcd_px, "1-based column index of x coordinate in the barcode file")
    LONG_INT_PARAM("icol-bcd-y", &in_icol_bcd_py, "1-based column index of y coordinate in the barcode file")
    LONG_INT_PARAM("icol-ftr-id", &in_icol_ftr_id, "1-based column index of feature ID in the barcode file")
    LONG_INT_PARAM("icol-ftr-name", &in_icol_ftr_name, "1-based column index of feature name in the barcode file")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_DOUBLE_PARAM("units-per-px", &units_per_px, "Coordinate unit per pixel")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (out_sgedir.empty() || in_sgedir.empty() || in_npy.empty())
        error("--in-npy, --in-sge and --out-sge must be specified");

    notice("Analysis started");

    // load the npy file
    try {
        npy::npy_data<int> npy_data = npy::read_npy<int>(in_npy);

        std::vector<unsigned long> shape = npy_data.shape;
        uint32_t nrow = shape[0];
        uint32_t ncol = shape[1];
        notice("Loaded npy file with %d rows and %d columns", nrow, ncol);

        // identify non-zero segments
        std::map<int32_t, std::vector<uint32_t>> segments;
        uint32_t nz = 0;
        for (uint32_t i = 0; i < nrow; i++) {
            for (uint32_t j = 0; j < ncol; j++) {
                int32_t seg = npy_data.data[i * ncol + j];
                if ( seg > 0 ) {
                    uint64_t key = ((uint64_t)i << 32 | (uint64_t)j);
                    segments[seg].push_back(key);
                    ++nz;
                }
            }
        }
        notice("Identified %zu segments with %d (%.5f) masked pixels", segments.size(), nz, nz / (double)nrow / ncol);        
    }
    catch (const std::exception& e) {
        error("Error reading npy file. Expecting a numpy array of integers: %s", e.what());
    }

    notice("Analysis finished");

    return 0;
}
