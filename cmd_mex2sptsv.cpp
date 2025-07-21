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

void write_sptsv_barcode(htsFile* wh_tsv, std::string& bcd, std::map<int32_t,int32_t>& iftr2cnts, std::vector<std::string>& ftr_combined, std::map<int32_t,int32_t>& iftr2oftr) {
    char buf[255];
    snprintf(buf, 255, "%08x", (uint32_t)rand() ); // 1-based index for the barcode
    int32_t sum_cnts = 0;
    for(std::map<int32_t,int32_t>::iterator it = iftr2cnts.begin(); it != iftr2cnts.end(); ++it) {
        int32_t iftr = it->first;
        int32_t cnt = it->second;
        sum_cnts += cnt;
    }
    hprintf(wh_tsv, "%s\t%s\t%zu\t%d", buf, bcd.c_str(), iftr2cnts.size(), sum_cnts);
    for(std::map<int32_t,int32_t>::iterator it = iftr2cnts.begin(); it != iftr2cnts.end(); ++it) {
        int32_t iftr = it->first;
        int32_t cnt = it->second;
        int32_t oftr = iftr2oftr[iftr];
        hprintf(wh_tsv, "\t%d %d", oftr, cnt);
    }
    hprintf(wh_tsv, "\n");
}


/////////////////////////////////////////////////////////////////////////
// dge2sptsv : Convert 10x Market Exchange (MEX) formatted DGE to Sparse TSV format in FICTURE2
////////////////////////////////////////////////////////////////////////
int32_t cmdMEX2SpTSV(int32_t argc, char **argv)
{
    std::string indir;
    std::string bcdf("barcodes.tsv.gz");
    std::string ftrf("features.tsv.gz");
    std::string mtxf("matrix.mtx.gz");
    int32_t icol_mtx = 3; // 1-based column index in the matrix file to use as the count
    int32_t seed = 0;  // random seed
    int32_t min_feature_count = 0;  // minimum feature count to include in the output
    std::string outprefix;
    std::string out_suffix_meta(".json");
    std::string out_suffix_tsv(".tsv");
    std::string out_suffix_feature_counts(".feature.counts.tsv");
    std::string colname_barcode_name("cell_id");
    std::string colname_random_key("random_key");
    std::string keyname_dictionary("dictionary");
    std::string keyname_header_info("header_info");
    std::string keyname_n_features("n_features");
    std::string keyname_n_modalities("n_modalities");
    std::string keyname_n_units("n_units");
    std::string keyname_offset_data("offset_data");
    bool use_gene_id = false;
    bool error_on_duplicate_gene_name = false;
    std::string include_ftr_list;
    std::string exclude_ftr_list;
    std::string include_bcd_list;
    std::string exclude_bcd_list;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input Options", NULL)
    LONG_STRING_PARAM("in-dir", &indir, "Input directory containing the MEX files")
    LONG_STRING_PARAM("bcd", &bcdf, "Barcode file name")
    LONG_STRING_PARAM("ftr", &ftrf, "Feature file name")
    LONG_STRING_PARAM("mtx", &mtxf, "Matrix file name")
    LONG_INT_PARAM("icol-mtx", &icol_mtx, "1-based column index in the matrix file to use as the count")

    LONG_PARAM_GROUP("Input Filtering Options", NULL)
    LONG_STRING_PARAM("include-feature-list", &include_ftr_list, "A file containing a list of input genes to be included (feature name of IDs)")
    LONG_STRING_PARAM("exclude-feature-list", &exclude_ftr_list, "A file containing a list of input genes to be excluded (feature name of IDs)")
    LONG_STRING_PARAM("include-barcode-list", &include_bcd_list, "A file containing a list of input barcode IDs to be included (feature name of IDs)")
    LONG_STRING_PARAM("exclude-barcode-list", &exclude_bcd_list, "A file containing a list of input barcode IDs to be excluded (feature name of IDs)")
    LONG_INT_PARAM("min-feature-count", &min_feature_count, "Minimum feature count to include in the output. Requires reading the input matrix twice")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_STRING_PARAM("out", &outprefix, "Output prefix")
    LONG_STRING_PARAM("out-suffix-meta", &out_suffix_meta, "Suffix for the metadata output file")
    LONG_STRING_PARAM("out-suffix-tsv", &out_suffix_tsv, "Suffix for the transcript output file")
    LONG_STRING_PARAM("out-suffix-feature-counts", &out_suffix_feature_counts, "Suffix for the feature counts output file")
    LONG_STRING_PARAM("colname-barcode-name", &colname_barcode_name, "Column name for barcode name")
    LONG_STRING_PARAM("colname-random-key", &colname_random_key, "Column name for random key")
    LONG_STRING_PARAM("keyname-dictionary", &keyname_dictionary, "Key name for the dictionary in the metadata file")
    LONG_STRING_PARAM("keyname-header-info", &keyname_header_info, "Key name for the header in the metadata file")
    LONG_STRING_PARAM("keyname-n-features", &keyname_n_features, "Key name for the number of features in the metadata file")
    LONG_STRING_PARAM("keyname-n-modalities", &keyname_n_modalities, "Key name for the number of modalities in the metadata file")
    LONG_STRING_PARAM("keyname-n-units", &keyname_n_units, "Key name for the number of units in the metadata file")
    LONG_STRING_PARAM("keyname-offset-data", &keyname_offset_data, "Key name for the offset data in the metadata file")

    LONG_PARAM_GROUP("Auxilary Output Options", NULL)
    LONG_INT_PARAM("seed", &seed, "Random seed for the random key generation")
    LONG_PARAM("use-gene-id", &use_gene_id, "Use gene ID instead of gene name in the output")
    LONG_PARAM("error-on-duplicate-gene-name", &error_on_duplicate_gene_name, "Report error on duplicate gene name in the feature file. If not set, the duplicate gene names will be renamed by combining with the gene IDs")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    // set random seed
    if ( seed > 0 ) {
        srand(seed);
    }
    else {
        srand(time(NULL));
    }

    if ( outprefix.empty() ) {
        error("--out must be specified");
    }

    if ( ! indir.empty() ) {
        bcdf = indir + "/" + bcdf;
        ftrf = indir + "/" + ftrf;
        mtxf = indir + "/" + mtxf;
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
    if ( stat(bcdf.c_str(), &sb) != 0 )
        error("Cannot find the barcode file %s", bcdf.c_str());
    if ( stat(ftrf.c_str(), &sb) != 0 )
        error("Cannot find the feature file %s", ftrf.c_str());
    if ( stat(mtxf.c_str(), &sb) != 0 )
        error("Cannot find the matrix file %s", mtxf.c_str());

    std::string out_tsv = outprefix + out_suffix_tsv;
    std::string out_meta = outprefix + out_suffix_meta;
    std::string out_feature_counts = outprefix + out_suffix_feature_counts;

    notice("Reading feature files...");
    // read the feature files
    std::vector<std::string> ftr_ids;
    std::vector<std::string> ftr_names;
    std::vector<std::string> ftr_combined;
    {
        tsv_reader ftr_tr(ftrf.c_str());
        std::set<std::string> ftr_id_set;
        std::map<std::string, int32_t> ftr_name2cnt;
        while (ftr_tr.read_line() > 0)
        {
            ftr_ids.push_back(ftr_tr.str_field_at(0));   // gene ID
            ftr_names.push_back(ftr_tr.str_field_at(1)); // gene name
            if ( !ftr_id_set.insert(ftr_tr.str_field_at(0)).second ) 
            {
                error("Duplicate gene ID found: %s", ftr_tr.str_field_at(0));
            }
            int32_t cnt = ++ftr_name2cnt[ftr_tr.str_field_at(1)];
            if ( error_on_duplicate_gene_name && cnt > 1 )
            {
                error("Duplicate gene name found: %s", ftr_tr.str_field_at(1));
            }
        }
        for(int32_t i = 0; i < (int32_t)ftr_ids.size(); ++i)
        {
            if ( use_gene_id )
                ftr_combined.push_back(ftr_ids[i]);
            else {
                if ( ftr_name2cnt[ftr_names[i]] > 1 ) {
                    ftr_combined.push_back(ftr_names[i] + "_" + ftr_ids[i]);
                }
                else {
                    ftr_combined.push_back(ftr_names[i]);
                }
            }
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
    }
    if (!exclude_bcd_list.empty())
    {
        tsv_reader bcd_tr(exclude_bcd_list.c_str());
        while (bcd_tr.read_line() > 0)
            exclude_bcd_set.insert(bcd_tr.str_field_at(0));
    }

    // load the exclude and include gene lists
    std::set<std::string> include_ftr_set;
    std::set<std::string> exclude_ftr_set;
    if (!include_ftr_list.empty())
    {
        tsv_reader ftr_tr(include_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
            include_ftr_set.insert(ftr_tr.str_field_at(0));
    }
    if (!exclude_ftr_list.empty())
    {
        tsv_reader ftr_tr(exclude_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
            exclude_ftr_set.insert(ftr_tr.str_field_at(0));
    }

    std::map<int32_t, int32_t> ftr2cnt; // total count for each feature
    if ( min_feature_count > 0 ) {
        // read the matrix file to get the total count for each feature
        notice("Reading the matrix file to get the total count for each feature...");
        tsv_reader mtx_tr(mtxf.c_str());
        while ( mtx_tr.read_line() ) {
            if ( mtx_tr.str_field_at(0)[0] == '%' ) { // skip the header lines
                continue;
            }
            int32_t iftr = mtx_tr.int_field_at(0); // feature index (1-based)
            int32_t cnt = mtx_tr.int_field_at(icol_mtx - 1); // count value (1-based)
            if ( cnt > 0 ) {
                ftr2cnt[iftr-1] += cnt; // accumulate the count for the feature
            }
        }
        mtx_tr.close();
    }

    // construct feature maps
    std::map<int32_t, int32_t> iftr2oftr; // index to map input feature index to output feature index (0-based)
    std::vector<int32_t> iftrs;
    for(int32_t i = 0; i < (int32_t)ftr_ids.size(); ++i)
    {
        if ( min_feature_count > 0 && (ftr2cnt[i] < min_feature_count ) ) {
            continue; // skip the feature if its total count is less than the minimum feature count
        }
        if ( include_ftr_set.empty() && exclude_ftr_set.empty() ) { // no filtering
            iftr2oftr[i] = (int32_t)iftrs.size(); // map input feature index to output feature index
            iftrs.push_back(i); // add to the output feature index
        }
        else if ( !include_ftr_set.empty() ) { // include list was specified
            if ( include_ftr_set.find(ftr_ids[i]) != include_ftr_set.end() ||
                 include_ftr_set.find(ftr_names[i]) != include_ftr_set.end() ) {
                iftr2oftr[i] = (int32_t)iftrs.size(); // add to the output feature index
                iftrs.push_back(i);
            }
        }
        else if ( !exclude_ftr_set.empty() ) { // exclude list was specified
            if ( exclude_ftr_set.find(ftr_ids[i]) == exclude_ftr_set.end() &&
                 exclude_ftr_set.find(ftr_names[i]) == exclude_ftr_set.end() ) {
                iftr2oftr[i] = (int32_t)iftrs.size(); // add to the output feature index
                iftrs.push_back(i);
            }
        }
        else {
            error("Something went wrong with the feature filtering options");
        }
    }

    tsv_reader bcd_tr(bcdf.c_str());
    std::vector<std::string> bcds;
    while( bcd_tr.read_line() > 0 )
    {
        bcds.push_back(bcd_tr.str_field_at(0));
    }
    bcd_tr.close();

    notice("Wrting output TSV file... %s", out_tsv.c_str());

    tsv_reader mtx_tr(mtxf.c_str());
    htsFile* wh_tsv = hts_open(out_tsv.c_str(), out_tsv.compare(out_tsv.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    mtx_tr.delimiter = ' ';
    int32_t nftrs = 0;
    int32_t nbcds = 0;
    uint64_t nlines = 0;
    int32_t cur_ibcd = -1;
    std::map<int32_t, int32_t> iftr2cnts;
    int32_t n_units = 0;
    while ( mtx_tr.read_line() ) {
        if ( mtx_tr.str_field_at(0)[0] == '%' ) { // skip the header lines
            continue;
        }
        if ( nlines == 0 ) {
            nftrs = mtx_tr.uint64_field_at(0); // number of features
            nbcds = mtx_tr.uint64_field_at(1); // number of barcodes
            nlines = mtx_tr.uint64_field_at(2); // number of lines
        }
        else {
            int32_t iftr = mtx_tr.int_field_at(0); // feature index (1-based)
            int32_t ibcd = mtx_tr.int_field_at(1); // barcode index (1-based)
            int32_t cnt = mtx_tr.int_field_at(icol_mtx - 1); // count value (1-based)
            if ( cnt > 0 && min_feature_count == 0 ) {
                ftr2cnt[iftr-1] += cnt;
            }
            if ( ibcd > cur_ibcd ) { // print the current barcode
                if ( ( cur_ibcd >= 0 ) && ( iftr2cnts.size() > 0 ) ) {
                    // check if the barcode is filtered or not
                    bool bcd_skip = ( ( ( !include_bcd_list.empty() ) && 
                                      ( include_bcd_set.find(bcds[ibcd-1]) == include_bcd_set.end() ) ) ||
                                      ( !exclude_bcd_list.empty() &&
                                       ( exclude_bcd_set.find(bcds[ibcd-1]) != exclude_bcd_set.end() ) ) );
                    if ( !bcd_skip) {
                        write_sptsv_barcode(wh_tsv, bcds[ibcd-1], iftr2cnts, ftr_combined, iftr2oftr);  
                        ++n_units;                  
                    }
                }
                iftr2cnts.clear();
                cur_ibcd = ibcd;
            }
            else if ( ibcd < cur_ibcd ) {
                error("Barcode is not stored in increasing order: %d < %d", ibcd, cur_ibcd);
                ++n_units;
            }
            if ( cnt > 0 ) {
                if ( iftr2oftr.find(iftr-1) != iftr2oftr.end() ) {
                    iftr2cnts[iftr-1] += cnt;
                }
            }
        }
    }
    mtx_tr.close();
    if ( iftr2cnts.size() > 0 ) {
        bool bcd_skip = ( ( ( !include_bcd_list.empty() ) && 
                    ( include_bcd_set.find(bcds[cur_ibcd-1]) == include_bcd_set.end() ) ) ||
                    ( !exclude_bcd_list.empty() &&
                    ( exclude_bcd_set.find(bcds[cur_ibcd-1]) != exclude_bcd_set.end() ) ) );
        if ( !bcd_skip ) {
            write_sptsv_barcode(wh_tsv, bcds[cur_ibcd-1], iftr2cnts, ftr_combined, iftr2oftr);
            ++n_units; 
        }
    }
    hts_close(wh_tsv);

    // write the feature counts
    notice("Writing feature counts...");
    htsFile* wh_counts = hts_open(out_feature_counts.c_str(), out_feature_counts.compare(out_feature_counts.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    for(int32_t i = 0; i < (int32_t)ftr_combined.size(); ++i) {
        hprintf(wh_counts, "%s\t%d\n", ftr_combined[i].c_str(), ftr2cnt[i]);
    }
    hts_close(wh_counts);

    notice("Writing metadata file... %s", out_meta.c_str());
    // write the metadata file
    htsFile* wh_meta = hts_open(out_meta.c_str(), out_meta.compare(out_meta.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    nlohmann::json j_meta;
    std::map<std::string, int32_t> j_dict;
    for(int32_t i=0; i < (int32_t)iftrs.size(); ++i) {
        j_dict[ftr_combined[iftrs[i]]] = iftr2oftr[iftrs[i]]; // 0-based index for the output feature
    }
    std::vector<std::string> j_header_info;
    j_header_info.push_back(colname_random_key);
    j_header_info.push_back(colname_barcode_name);

    j_meta[keyname_dictionary] = j_dict;
    j_meta[keyname_header_info] = j_header_info;
    j_meta[keyname_n_features] = (int32_t)iftrs.size(); // number of features
    j_meta[keyname_n_modalities] = 1; // number of modalities
    j_meta[keyname_offset_data] = 2; // offset for the data columns
    j_meta[keyname_n_units] = n_units; // number of units
    hprintf(wh_meta, "%s\n", j_meta.dump(4).c_str());
    hts_close(wh_meta);

    notice("Analysis finished");

    return 0;
}
