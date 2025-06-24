#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include "nlohmann/json.hpp"

/////////////////////////////////////////////////////////////////////////
// mex_subset : Subset 10x Market Exchange (MEX) formatted to specific 
////////////////////////////////////////////////////////////////////////
int32_t cmdMEXSubset(int32_t argc, char **argv)
{
    std::string in_dir;
    std::string in_bcdf("barcodes.tsv.gz");
    std::string in_ftrf("features.tsv.gz");
    std::string in_mtxf("matrix.mtx.gz");
    int32_t icol_mtx = 3; // 1-based column index in the matrix file to use as the count
    int32_t min_feature_count = 0;  // minimum feature count to include in the output
    int32_t min_barcode_count = 0;  // minimum barcode count to include in the output
    std::string out_dir;
    std::string out_bcdf("barcodes.tsv.gz");
    std::string out_ftrf("features.tsv.gz");
    std::string out_mtxf("matrix.mtx.gz");
    std::string include_ftr_list;
    std::string exclude_ftr_list;
    std::string include_bcd_list;
    std::string exclude_bcd_list;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input Options", NULL)
    LONG_STRING_PARAM("in-dir", &in_dir, "Input directory containing the MEX files")
    LONG_STRING_PARAM("in-bcd", &in_bcdf, "Input Barcode file name")
    LONG_STRING_PARAM("in-ftr", &in_ftrf, "Input Feature file name")
    LONG_STRING_PARAM("in-mtx", &in_mtxf, "Input Matrix file name")
    LONG_INT_PARAM("icol-mtx", &icol_mtx, "1-based column index in the matrix file to use as the count")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_STRING_PARAM("out-dir", &out_dir, "Output directory containing the MEX files")
    LONG_STRING_PARAM("out-bcd", &out_bcdf, "Output Barcode file name")
    LONG_STRING_PARAM("out-ftr", &out_ftrf, "Output Feature file name")
    LONG_STRING_PARAM("out-mtx", &out_mtxf, "Output Matrix file name")

    LONG_PARAM_GROUP("Input Filtering Options", NULL)
    LONG_STRING_PARAM("include-feature-list", &include_ftr_list, "A file containing a list of input genes to be included (feature name of IDs)")
    LONG_STRING_PARAM("exclude-feature-list", &exclude_ftr_list, "A file containing a list of input genes to be excluded (feature name of IDs)")
    LONG_STRING_PARAM("include-barcode-list", &include_bcd_list, "A file containing a list of input barcode IDs to be included (feature name of IDs)")
    LONG_STRING_PARAM("exclude-barcode-list", &exclude_bcd_list, "A file containing a list of input barcode IDs to be excluded (feature name of IDs)")
    LONG_INT_PARAM("min-feature-count", &min_feature_count, "Minimum feature count to include in the output. Requires reading the input matrix twice")
    LONG_INT_PARAM("min-barcode-count", &min_barcode_count, "Minimum barcode count to include in the output. Requires reading the input matrix twice")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if ( ! in_dir.empty() ) {
        in_bcdf = in_dir + "/" + in_bcdf;
        in_ftrf = in_dir + "/" + in_ftrf;
        in_mtxf = in_dir + "/" + in_mtxf;
    }

    if ( !out_dir.empty() ) {
        out_bcdf = out_dir + "/" + out_bcdf;
        out_ftrf = out_dir + "/" + out_ftrf;
        out_mtxf = out_dir + "/" + out_mtxf;
        if (makePath(out_dir))
        {
            notice("Successfully created the directory %s", out_dir.c_str());
        }
        else
        {
            notice("The directory %s already exists", out_dir.c_str());
        }
    }

    if ( ( !include_ftr_list.empty() ) && ( !exclude_ftr_list.empty() ) )
    {
        error("Cannot specify both --include-feature-list and --exclude-feature-list");
    }
    if ( ( !include_bcd_list.empty() ) && ( !exclude_bcd_list.empty() ) )
    {
        error("Cannot specify both --include-barcode-list and --exclude-barcode-list");
    }

    // check if the input files exist
    struct stat sb;
    if ( stat(in_bcdf.c_str(), &sb) != 0 )
        error("Cannot find the barcode file %s", in_bcdf.c_str());
    if ( stat(in_ftrf.c_str(), &sb) != 0 )
        error("Cannot find the feature file %s", in_ftrf.c_str());
    if ( stat(in_mtxf.c_str(), &sb) != 0 )
        error("Cannot find the matrix file %s", in_mtxf.c_str());

    notice("Reading feature files...");
    // read the feature files
    std::vector<std::string> ftr_ids;
    std::vector<std::string> ftr_names;
    std::vector<std::string> ftr_types; // feature types, e.g. "Gene Expression"
    std::map<std::string, std::vector<int32_t> > ftr_name2idx;
    std::set<std::string> ftr_id_set; // to check for duplicate feature IDs
    {
        tsv_reader ftr_tr(in_ftrf.c_str());
        ftr_tr.delimiter = '\t'; // set the delimiter to tab
        std::set<std::string> ftr_id_set;
        while (ftr_tr.read_line())
        {
            const char* id = ftr_tr.str_field_at(0);
            const char* name = ftr_tr.str_field_at(1);
            ftr_ids.push_back(id);   // gene ID
            ftr_names.push_back(name); // gene name
            if ( ftr_tr.nfields >= 3 ) {
                const char* type = ftr_tr.str_field_at(2);
                ftr_types.push_back(type); // gene type
            }
            else {
                ftr_types.push_back("Gene Expression"); // default type if not specified
            }
            if ( ftr_id_set.find(id) != ftr_id_set.end() ) {
                error("Duplicate feature ID found: %s", id);
            }
            ftr_id_set.insert(id); // add the feature ID to the set
            ftr_name2idx[name].push_back((int32_t)ftr_ids.size() - 1); // store the index of the gene ID
        }
        ftr_tr.close();
    }

    notice("Processing barcode and feature filtering options....");
    // load the exclude and include barcode lists
    std::set<std::string> include_bcd_set;
    std::set<std::string> exclude_bcd_set;
    if (!include_bcd_list.empty())
    {
        tsv_reader bcd_tr(include_bcd_list.c_str());
        while (bcd_tr.read_line() > 0)
            include_bcd_set.insert(bcd_tr.str_field_at(0));
        bcd_tr.close();
        notice("Loaded %zu barcodes from the include list", include_bcd_set.size());
    }
    if (!exclude_bcd_list.empty())
    {
        tsv_reader bcd_tr(exclude_bcd_list.c_str());
        while (bcd_tr.read_line() > 0)
            exclude_bcd_set.insert(bcd_tr.str_field_at(0));
        bcd_tr.close();
        notice("Loaded %zu barcodes from the exclude list", exclude_bcd_set.size());
    }

    // load the exclude and include gene lists
    std::set<std::string> include_ftr_set;
    std::set<std::string> exclude_ftr_set;
    if (!include_ftr_list.empty())
    {
        tsv_reader ftr_tr(include_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
            include_ftr_set.insert(ftr_tr.str_field_at(0));
        ftr_tr.close();
        notice("Loaded %zu features from the include list", include_ftr_set.size());
    }
    if (!exclude_ftr_list.empty())
    {
        tsv_reader ftr_tr(exclude_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
            exclude_ftr_set.insert(ftr_tr.str_field_at(0));
        ftr_tr.close();
        notice("Loaded %zu features from the exclude list", exclude_ftr_set.size());
    }

    std::map<int32_t, int32_t> ftr2cnt; // total count for each feature
    std::map<int32_t, int32_t> bcd2cnt; // total count for each barcode
    if ( min_feature_count > 0 || min_barcode_count > 0 ) 
    {
        // read the matrix file to get the total count for each feature
        notice("Reading the matrix file to get the total count for each feature...");
        tsv_reader mtx_tr(in_mtxf.c_str());
        uint64_t nlines = 0;
        uint64_t nnz = 0;
        while ( mtx_tr.read_line() ) {
            if ( mtx_tr.str_field_at(0)[0] == '%' ) { // skip the header lines
                continue;
            }
            if ( nlines == 0 ) {
                nnz = mtx_tr.uint64_field_at(2); // total number of non-zero entries
                ++nlines;
            }
            else {
                int32_t iftr = mtx_tr.int_field_at(0); // feature index (1-based)
                int32_t ibcd = mtx_tr.int_field_at(1); // barcode index (1-based)
                int32_t cnt = mtx_tr.int_field_at(icol_mtx - 1); // count value (1-based)
                if ( cnt > 0 ) {
                    ftr2cnt[iftr-1] += cnt; // accumulate the count for the feature
                    bcd2cnt[ibcd-1] += cnt; // accumulate the count for the barcode
                }
                ++nlines;
                if ( nlines % 1000000 == 0 ) {
                    notice("Processed %llu lines, total non-zero entries: %llu (%.5lf)", nlines, nnz, nlines / (double)nnz);
                }
            }
        }
        mtx_tr.close();
    }

    // construct feature maps
    std::map<int32_t, int32_t> iftr2oftr; // index to map input feature index to output feature index (0-based)
    std::vector<int32_t> iftrs;
    htsFile* wh_ftr = hts_open(out_ftrf.c_str(), out_ftrf.compare(out_ftrf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wh_ftr == NULL ) {
        error("Cannot open feature file %s for writing", out_ftrf.c_str());
    }
    for(int32_t i = 0; i < (int32_t)ftr_ids.size(); ++i)
    {
        if ( min_feature_count > 0 && (ftr2cnt[i] < min_feature_count ) ) {
            continue; // skip the feature if its total count is less than the minimum feature count
        }
        if ( include_ftr_set.empty() && exclude_ftr_set.empty() ) { // no filtering
            iftr2oftr[i] = (int32_t)iftrs.size(); // map input feature index to output feature index
            iftrs.push_back(i); // add to the output feature index
            hprintf(wh_ftr, "%s\t%s\t%s\n", ftr_ids[i].c_str(), ftr_names[i].c_str(), ftr_types[i].c_str()); // write feature index (1-based) and name
        }
        else if ( !include_ftr_set.empty() ) { // include list was specified
            if ( include_ftr_set.find(ftr_ids[i]) != include_ftr_set.end() ||
                 include_ftr_set.find(ftr_names[i]) != include_ftr_set.end() ) {
                iftr2oftr[i] = (int32_t)iftrs.size(); // add to the output feature index
                iftrs.push_back(i);
                hprintf(wh_ftr, "%s\t%s\t%s\n", ftr_ids[i].c_str(), ftr_names[i].c_str(), ftr_types[i].c_str()); // write feature index (1-based) and name
            }
        }
        else if ( !exclude_ftr_set.empty() ) { // exclude list was specified
            if ( exclude_ftr_set.find(ftr_ids[i]) == exclude_ftr_set.end() &&
                 exclude_ftr_set.find(ftr_names[i]) == exclude_ftr_set.end() ) {
                iftr2oftr[i] = (int32_t)iftrs.size(); // add to the output feature index
                iftrs.push_back(i);
                hprintf(wh_ftr, "%s\t%s\t%s\n", ftr_ids[i].c_str(), ftr_names[i].c_str(), ftr_types[i].c_str()); // write feature index (1-based) and name
            }
        }
        else {
            error("Something went wrong with the feature filtering options");
        }
    }
    hts_close(wh_ftr);

    tsv_reader bcd_tr(in_bcdf.c_str());
    std::vector<std::string> bcds;
    while( bcd_tr.read_line() > 0 )
    {
        bcds.push_back(bcd_tr.str_field_at(0));
    }
    bcd_tr.close();

    std::set<int32_t> iftrs_set;
    std::set<int32_t> ibcds_set;
    uint64_t nlines = 0;
    {
        // read the matrix file to get the total count for each feature
        notice("Reading the matrix file for the 2nd time to get the total number of lines..");
        tsv_reader mtx_tr(in_mtxf.c_str());
        mtx_tr.delimiter = ' ';
        uint64_t nnz = 0;
        uint64_t nproc = 0;
        while ( mtx_tr.read_line() ) {
            if ( mtx_tr.str_field_at(0)[0] == '%' ) { // skip the header lines
                continue;
            }
            if ( nproc == 0 ) {
                nnz = mtx_tr.uint64_field_at(2); // total number of non-zero entries
                ++nproc;
            }
            else {
                int32_t iftr = mtx_tr.int_field_at(0); // feature index (1-based)
                int32_t ibcd = mtx_tr.int_field_at(1); // barcode index (1-based)
                ++nproc;
                if ( nproc % 1000000 == 0 ) {
                    notice("Processed %llu lines, passed %llu, total non-zero entries: %llu (%.5lf)", nproc, nlines, nnz, nproc / (double)nnz);
                }
                if ( iftr2oftr.find(iftr-1) == iftr2oftr.end() ) {
                    continue; // skip the feature if it is not in the output feature list
                }
                if ( min_feature_count > 0 && ftr2cnt[iftr-1] < min_feature_count ) {
                    continue; // skip the feature if its total count is less than the minimum feature count
                }
                if ( min_barcode_count > 0 && bcd2cnt[ibcd-1] < min_barcode_count ) {
                    continue; // skip the barcode if its total count is less than the minimum barcode count
                }
                if ( ( !include_bcd_list.empty() ) && 
                        ( include_bcd_set.find(bcds[ibcd-1]) == include_bcd_set.end() ) ) {
                        continue; // skip the barcode if it is not in the include list
                }
                else if ( !exclude_bcd_list.empty() &&
                        ( exclude_bcd_set.find(bcds[ibcd-1]) != exclude_bcd_set.end() ) ) {
                    continue; // skip the barcode if it is in the exclude list
                }
                ++nlines; // count the number of lines
                iftrs_set.insert(iftr-1); // add the feature index to the set
                ibcds_set.insert(ibcd-1); // add the barcode index to the set
            }
        }
        mtx_tr.close();
    }

    htsFile* wh_mtx = hts_open(out_mtxf.c_str(), out_mtxf.compare(out_mtxf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wh_mtx == NULL ) {
        error("Cannot open matrix file %s for writing", out_mtxf.c_str());
    }
    htsFile* wh_bcd = hts_open(out_bcdf.c_str(), out_bcdf.compare(out_bcdf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wh_bcd == NULL ) {
        error("Cannot open barcode file %s for writing", out_bcdf.c_str());
    }
    uint64_t nlines2 = 0;
    {
        tsv_reader mtx_tr(in_mtxf.c_str());
        int32_t prev_ibcd = -1;
        uint64_t nnz = 0, nproc = 0;
        notice("Reading the matrix file for the 2rd time to write the actual output..");
        while ( mtx_tr.read_line() ) {
            if ( mtx_tr.str_field_at(0)[0] == '%' ) { // skip the header lines
                for(int32_t i = 0; i < mtx_tr.nfields; ++i) {
                    if ( i > 0 ) hprintf(wh_mtx, " ");
                    hprintf(wh_mtx, "%s", mtx_tr.str_field_at(i));
                }
                hprintf(wh_mtx, "\n");
                continue;
            }
            if ( nproc == 0 ) {            
                hprintf(wh_mtx, "%zu %zu %zu\n", iftrs_set.size(), ibcds_set.size(), nlines-1);
                ++nproc;
                nnz = mtx_tr.uint64_field_at(2); // total number of non-zero entries
            }
            else {
                int32_t iftr = mtx_tr.int_field_at(0); // feature index (1-based)
                int32_t ibcd = mtx_tr.int_field_at(1); // barcode index (1-based)
                int32_t cnt = mtx_tr.int_field_at(icol_mtx - 1); // count value (1-based)
                ++nproc;
                if ( nlines2 % 1000000 == 0 ) {
                    notice("Processed %llu lines, passed %llu, total non-zero entries: %llu (%.5lf)", nproc, nlines2, nnz, nproc / (double)nnz);
                }
                if ( iftr2oftr.find(iftr-1) == iftr2oftr.end() ) {
                    continue; // skip the feature if it is not in the output feature list
                }
                if ( min_feature_count > 0 && ftr2cnt[iftr-1] < min_feature_count ) {
                    continue; // skip the feature if its total count is less than the minimum feature count
                }
                if ( min_barcode_count > 0 && bcd2cnt[ibcd-1] < min_barcode_count ) {
                    continue; // skip the barcode if its total count is less than the minimum barcode count
                }
                if ( ( !include_bcd_list.empty() ) && 
                        ( include_bcd_set.find(bcds[ibcd-1]) == include_bcd_set.end() ) ) {
                        continue; // skip the barcode if it is not in the include list
                }
                else if ( !exclude_bcd_list.empty() &&
                        ( exclude_bcd_set.find(bcds[ibcd-1]) != exclude_bcd_set.end() ) ) {
                    continue; // skip the barcode if it is in the exclude list
                }
                if ( ibcd != prev_ibcd ) {
                    hprintf(wh_bcd, "%s\n", bcds[ibcd-1].c_str()); // write the barcode 
                    prev_ibcd = ibcd; // update the previous barcode index
                }
                hprintf(wh_mtx, "%d %d %d\n", iftr2oftr[iftr-1] + 1, ibcd, cnt); // write the feature index (1-based), barcode index (1-based) and count
                ++nlines2; // count the number of lines
            }
        }
        mtx_tr.close();
    }
    hts_close(wh_mtx);
    hts_close(wh_bcd);


    notice("Analysis finished");

    return 0;
}
