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
#include <fstream>
#include "nlohmann/json.hpp"


/////////////////////////////////////////////////////////////////////////
// sptsv2mex : Convert Sparse TSV format in FICTURE2 to 10x Market Exchange (MEX) formatted DGE 
////////////////////////////////////////////////////////////////////////
int32_t cmdSpTSV2MEX(int32_t argc, char **argv)
{
    std::string intsv;
    std::string inmeta;
    std::string outdir; // output directory
    std::string bcdf("barcodes.tsv.gz");
    std::string ftrf("features.tsv.gz");
    std::string mtxf("matrix.mtx.gz");
    std::string colname_random_key("random_key");
    bool include_random_key = false;
    std::string keyname_dictionary("dictionary");
    std::string keyname_header_info("header_info");
    std::string keyname_n_features("n_features");
    std::string keyname_n_units("n_units");
    std::string keyname_offset_data("offset_data");
    std::string delim(":");
    std::string feature_type("Gene Expression"); // default feature type

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input Options", NULL)
    LONG_STRING_PARAM("tsv", &intsv, "Input TSV file containing the sparse matrix data")
    LONG_STRING_PARAM("json", &inmeta, "Input JSON metadata file containing the header information")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_STRING_PARAM("out-dir", &outdir, "Output directory for the MEX files")
    LONG_STRING_PARAM("bcd", &bcdf, "Barcode file name")
    LONG_STRING_PARAM("ftr", &ftrf, "Feature file name")
    LONG_STRING_PARAM("mtx", &mtxf, "Matrix file name")

    LONG_PARAM_GROUP("Auxiliary Input/Output Options", NULL)
    LONG_STRING_PARAM("colname-random-key", &colname_random_key, "Column name for the random key in the output")
    LONG_PARAM("include-random-key", &include_random_key, "Include the random key in the output")
    LONG_STRING_PARAM("keyname-dictionary", &keyname_dictionary, "Key name for the dictionary in the metadata file")
    LONG_STRING_PARAM("keyname-header-info", &keyname_header_info, "Key name for the header information in the metadata file")
    LONG_STRING_PARAM("keyname-n-features", &keyname_n_features, "Key name for the number of features in the metadata file")
    LONG_STRING_PARAM("keyname-n-units", &keyname_n_units, "Key name for the number of units in the metadata file")
    LONG_STRING_PARAM("keyname-offset-data", &keyname_offset_data, "Key name for the offset data in the metadata file")
    LONG_STRING_PARAM("delim", &delim, "Delimiter to use as barcode names when multiple headers fields are combined")
    LONG_STRING_PARAM("feature-type", &feature_type, "Feature type to use in the output files (default: Gene Expression)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if ( outdir.empty() ) {
        error("--out must be specified");
    }

    if ( ! outdir.empty() ) {
        bcdf = outdir + "/" + bcdf;
        ftrf = outdir + "/" + ftrf;
        mtxf = outdir + "/" + mtxf;
        if (makePath(outdir))
        {
            notice("Successfully created the directory %s", outdir.c_str());
        }
        else
        {
            notice("The directory %s already exists", outdir.c_str());
        }
    }

    notice("Reading metadata file %s", inmeta.c_str());
    nlohmann::json json_data;
    int32_t icol_random_key = -1; // column index for the random key, -1 if not present
    int32_t offset_data = -1;     // default offset for the data in the input TSV file
    int32_t n_features = -1; // number of features in the input TSV file
    int32_t n_units = -1;    // number of units in the input TSV file
    std::vector<std::string> feature_names;
    {
        std::ifstream meta_file(inmeta);
        if (!meta_file.is_open()) {
            error("Cannot open metadata file %s", inmeta.c_str());
        }
        meta_file >> json_data;

        try {
            std::vector<std::string> hdrs = json_data[keyname_header_info].get<std::vector<std::string> >();
            for(int32_t i=0; i < (int32_t)json_data[keyname_header_info].size(); ++i) {
                if ( hdrs[i].compare(colname_random_key) == 0 ) {
                    icol_random_key = i;
                }
            }
        }
        catch (nlohmann::json::exception& e) {
            error("Cannot read header information from the metadata file %s: %s", inmeta.c_str(), e.what());
        }
        offset_data = json_data["offset_data"].get<int32_t>();
        n_features = json_data[keyname_n_features].get<int32_t>();
        n_units = json_data[keyname_n_units].get<int32_t>();
        try {
            std::map<std::string, int32_t> dict = json_data[keyname_dictionary].get<std::map<std::string, int32_t> >();
            feature_names.resize(n_features);
            for (const auto& kv : dict) {
                if ( kv.second < 0 || kv.second >= n_features ) {
                    error("Invalid feature index %d for feature %s in the metadata file %s", kv.second, kv.first.c_str(), inmeta.c_str());
                }
                feature_names[kv.second] = kv.first; // map feature index to feature name
            }
            // check if all features are present
            for(int32_t i = 0; i < n_features; ++i) {
                if ( feature_names[i].empty() ) {
                    error("Feature %d is not present in the metadata file %s", i, inmeta.c_str());
                }
            }
        }
        catch (nlohmann::json::exception& e) {
            error("Cannot read dictionary from the metadata file %s: %s", inmeta.c_str(), e.what());
        }

        if ( icol_random_key < 0 && include_random_key ) {
            error("Random key column %s not found in the metadata file %s", colname_random_key.c_str(), inmeta.c_str());
        }
        if ( offset_data < 0 ) {
            error("Invalid offset_data value %d in the metadata file %s", offset_data, inmeta.c_str());
        }
        if ( n_features <= 0 ) {
            error("Invalid number of features %d in the metadata file %s", n_features, inmeta.c_str());
        }
        if ( n_units <= 0 ) {
            error("Invalid number of units %d in the metadata file %s", n_units, inmeta.c_str());
        }
    }

    htsFile* wf = hts_open(ftrf.c_str(), ftrf.compare(ftrf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open feature file %s for writing", ftrf.c_str());
    }
    for(int32_t i = 0; i < n_features; ++i)
    {
        hprintf(wf, "%s\t%s\t%s\n", feature_names[i].c_str(), feature_names[i].c_str(), feature_type.c_str()); // write feature index (1-based) and name
    }
    hts_close(wf);

    // computing the total number of lines by reading the input TSV file
    uint64_t nnz = 0;
    {
        tsv_reader tsv_tr;
        int32_t nlines = 0;
        if ( !tsv_tr.open(intsv.c_str()) ) {
            error("Cannot open input TSV file %s for reading", intsv.c_str());
        }
        while (tsv_tr.read_line()) {
            nnz += tsv_tr.int_field_at(offset_data);
            ++nlines;
        }
        tsv_tr.close();
        notice("Total number of non-zero entries in the input TSV file: %llu", nnz);
        if ( nlines != n_units ) {
            error("The number of lines in the input TSV file (%d) does not match the number of units (%d) in the metadata file", nlines, n_units);
        }
    }

    htsFile* wf_bcd = hts_open(bcdf.c_str(), bcdf.compare(bcdf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wf_bcd == NULL ) {
        error("Cannot open barcode file %s for writing", bcdf.c_str());
    }
    htsFile* wf_mtx = hts_open(mtxf.c_str(), mtxf.compare(mtxf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wf_mtx == NULL ) {
        error("Cannot open matrix file %s for writing", mtxf.c_str());
    }
    hprintf(wf_mtx, "%%%%MatrixMarket matrix coordinate real general\n");
    hprintf(wf_mtx, "%%\n");
    hprintf(wf_mtx, "%d %d %llu\n", n_features, n_units, nnz); // write the header for the matrix file

    tsv_reader tsv_tr;
    if ( !tsv_tr.open(intsv.c_str()) ) {
        error("Cannot open input TSV file %s for reading", intsv.c_str());
    }   
    int32_t ibcd = 0;
    while( tsv_tr.read_line() ) {
        std::string barcode;
        for(int32_t i = 0; i < offset_data; ++i) {
            if ( ( i == icol_random_key ) && ( !include_random_key ) ) {
                continue; // skip the random key if not included
            }
            if ( barcode.empty() ) {
                barcode.assign(tsv_tr.str_field_at(i));
            }
            else {
                barcode += (delim + tsv_tr.str_field_at(i)); // combine multiple fields with the specified delimiter
            }
        }
        ++ibcd;
        hprintf(wf_bcd, "%s\n", barcode.c_str()); // write the barcode to the barcode file
        int32_t n_genes = tsv_tr.int_field_at(offset_data);
        int32_t n_mols = tsv_tr.int_field_at(offset_data + 1);
        if ( tsv_tr.nfields != 2 * n_genes + offset_data + 2 ) {
            error("The number of fields in the input TSV file (%d) does not match the expected number (%d)", tsv_tr.nfields, n_genes + offset_data + 2);
        }
        int32_t sum_cnt = 0;
        for(int32_t i = 0; i < n_genes; ++i) {
            const char* s = tsv_tr.str_field_at(2*i + offset_data + 2);
            const char* t = tsv_tr.str_field_at(2*i + offset_data + 3);
            if ( t == NULL ) {
                error("Invalid data format in the input TSV file %s at line %d: expected space-separated values", intsv.c_str(), tsv_tr.nlines);
            }
            int32_t iftr = atoi(s);
            int32_t cnt = atoi(t);
            sum_cnt += cnt;
            hprintf(wf_mtx, "%d %d %d\n", iftr + 1, ibcd, cnt); // write the matrix entry (1-based index)
        }
        if ( sum_cnt != n_mols ) {
            error("The sum of counts (%d) does not match the number of molecules (%d) in the input TSV file %s at line %d", sum_cnt, n_mols, intsv.c_str(), tsv_tr.nlines);
        }
    }
    tsv_tr.close();
    if ( hts_close(wf_bcd) != 0 ) {
        error("Cannot close barcode file %s", bcdf.c_str());
    }
    if ( hts_close(wf_mtx) != 0 ) {
        error("Cannot close matrix file %s", mtxf.c_str());
    }

    notice("Analysis finished");

    return 0;
}
