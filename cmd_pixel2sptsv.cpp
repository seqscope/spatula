#include "nlohmann/json.hpp"
#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"
#include "spatula.h"
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>

/////////////////////////////////////////////////////////////////////////
// pixel2sptsv : Convert pixel-level transcript data into sparse TSV format
//               segmented by cells
////////////////////////////////////////////////////////////////////////
int32_t cmdPixel2SpTSV(int32_t argc, char **argv) {
    std::string pixelf;     // pixel-level transcript data
    std::string in_col_id("cell_id");  // column name in the input file for cell ID
    std::string ignore_ids("UNASSIGNED,-1,0,NA"); // IDs to be considered as null and ignored
    std::string in_col_ftr("gene");    // column name in the input file for feature name
    std::string in_col_cnt("count");   // column name in the input file for count  
    std::string in_col_x("X"); // column name in the input file for x coordinate
    std::string in_col_y("Y"); // column name in the input file for y coordinate
    int32_t idx_col_x = 1;   // 1-based index for x coordinate column
    int32_t idx_col_y = 2;   // 1-based index for y coordinate column
    int32_t idx_col_ftr = 3; // 1-based index for feature name column
    int32_t idx_col_cnt = 4; // 1-based index for count column
    int32_t idx_col_id = 5; // 1-based index for cell ID column

    bool skip_cell_tsv = false; // skip generating cell-level TSV file
    bool add_xy = false;        // add X and Y coordinate to the output TSV file
    bool no_header = false;     // if true, input file has no header line

    int32_t min_feature_count = 1; // minimum feature count to include in the output
    int32_t min_cell_count = 1;    // minimum cell count to include in the output, after filtering by features
    std::string outprefix; // prefix of sparse TSV output files 
    std::string out_suffix_json(".json");
    std::string out_suffix_tsv(".tsv");
    std::string out_suffix_cell_meta(".cell.metadata.tsv");
    std::string out_suffix_feature_counts(".feature.counts.tsv");
    std::string out_col_id("cell_id"); // column name in the output file for cell ID
    std::string out_col_random_key("random_key");
    std::string out_col_x("X"); // column name in the output file for x coordinate
    std::string out_col_y("Y"); // column name in the output file for y coordinate
    std::string keyname_dictionary("dictionary");
    std::string keyname_header_info("header_info");
    std::string keyname_n_features("n_features");
    std::string keyname_n_modalities("n_modalities");
    std::string keyname_n_units("n_units");
    std::string keyname_offset_data("offset_data");
    std::string include_ftr_list;
    std::string exclude_ftr_list;
    std::string include_id_list;
    std::string exclude_id_list;

    int32_t seed = 0;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input Options", NULL)
    LONG_STRING_PARAM("pixel", &pixelf, "Input pixel-level transcript data")
    LONG_STRING_PARAM("in-col-id", &in_col_id, "Column name in the input file for cell ID")
    LONG_STRING_PARAM("in-col-ftr", &in_col_ftr, "Column name in the input file for feature name")
    LONG_STRING_PARAM("in-col-cnt", &in_col_cnt, "Column name in the input file for count")
    LONG_STRING_PARAM("in-col-x", &in_col_x, "Column name in the input file for x coordinate")
    LONG_STRING_PARAM("in-col-y", &in_col_y, "Column name in the input file for y coordinate")
    LONG_INT_PARAM("idx-col-x", &idx_col_x, "1-based index for x coordinate column, effective only with --no-header")
    LONG_INT_PARAM("idx-col-y", &idx_col_y, "1-based index for y coordinate columnm, effective only with --no-header")
    LONG_INT_PARAM("idx-col-ftr", &idx_col_ftr, "1-based index for feature name columnm, effective only with --no-header")
    LONG_INT_PARAM("idx-col-cnt", &idx_col_cnt, "1-based index for count column, effective only with --no-header")
    LONG_INT_PARAM("idx-col-id", &idx_col_id, "1-based index for cell ID column, effective only with --no-header")
    LONG_PARAM("no-header", &no_header, "If set, the input file has no header line")

    LONG_PARAM_GROUP("Input Filtering Options", NULL)
    LONG_STRING_PARAM("ignore-ids", &ignore_ids, "IDs to be considered as null and ignored")
    LONG_STRING_PARAM("include-feature-list", &include_ftr_list, "A file containing a list of input genes to be included (feature name of IDs)")
    LONG_STRING_PARAM("exclude-feature-list", &exclude_ftr_list, "A file containing a list of input genes to be excluded (feature name of IDs)")
    LONG_STRING_PARAM("include-id-list", &include_id_list, "A file containing a list of input IDs to be included")
    LONG_STRING_PARAM("exclude-id-list", &exclude_id_list, "A file containing a list of input IDs to be excluded")
    LONG_INT_PARAM("min-cell-count", &min_cell_count, "Minimum cell count to include in the output.")
    LONG_INT_PARAM("min-feature-count", &min_feature_count, "Minimum feature count to include in the output.")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_STRING_PARAM("out", &outprefix, "Output prefix")
    LONG_STRING_PARAM("out-suffix-json", &out_suffix_json, "Suffix for the metadata JSON file")
    LONG_STRING_PARAM("out-suffix-cell-meta", &out_suffix_cell_meta, "Suffix for the cell metadata output file")
    LONG_STRING_PARAM("out-suffix-tsv", &out_suffix_tsv, "Suffix for the transcript output file")
    LONG_STRING_PARAM("out-suffix-feature-counts", &out_suffix_feature_counts, "Suffix for the feature counts output file")
    LONG_STRING_PARAM("out-col-id", &out_col_id, "Column name for cell ID")
    LONG_STRING_PARAM("out-col-random-key", &out_col_random_key, "Column name for random key")
    LONG_STRING_PARAM("out-col-x", &out_col_x, "Column name for x coordinate")
    LONG_STRING_PARAM("out-col-y", &out_col_y, "Column name for y coordinate")
    LONG_STRING_PARAM("keyname-dictionary", &keyname_dictionary, "Key name for the dictionary in the metadata file")
    LONG_STRING_PARAM("keyname-header-info", &keyname_header_info, "Key name for the header in the metadata file")
    LONG_STRING_PARAM("keyname-n-features", &keyname_n_features, "Key name for the number of features in the metadata file")
    LONG_STRING_PARAM("keyname-n-modalities", &keyname_n_modalities, "Key name for the number of modalities in the metadata file")
    LONG_STRING_PARAM("keyname-n-units", &keyname_n_units, "Key name for the number of units in the metadata file")
    LONG_STRING_PARAM("keyname-offset-data", &keyname_offset_data, "Key name for the offset data in the metadata file")
    LONG_PARAM("skip-cell-tsv", &skip_cell_tsv, "Skip generating cell TSV output")
    LONG_PARAM("add-xy", &add_xy, "Add x and y coordinates to the output sparse TSV file")

    LONG_PARAM_GROUP("Auxilary Output Options", NULL)
    LONG_INT_PARAM("seed", &seed, "Random seed for the random key generation")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    // set random seed
    if (seed > 0) {
        srand(seed);
    } else {
        srand(time(NULL));
    }

    if (outprefix.empty()) {
        error("--out must be specified");
    }

    if ((!include_ftr_list.empty()) && (!exclude_ftr_list.empty())) {
        error("Cannot specify both --include-feature-list and --exclude-feature-list");
    }
    if ((!include_id_list.empty()) && (!exclude_id_list.empty())) {
        error("Cannot specify both --include-id-list and --exclude-id-list");
    }

    // check if the input files exist
    struct stat sb;
    if (stat(pixelf.c_str(), &sb) != 0)
        error("Cannot find the pixel file %s", pixelf.c_str());

    std::string out_tsv = outprefix + out_suffix_tsv;
    std::string out_json = outprefix + out_suffix_json;
    std::string out_feature_counts = outprefix + out_suffix_feature_counts;
    std::string out_cell_meta = outprefix + out_suffix_cell_meta;
    bool compute_xy = (add_xy || !skip_cell_tsv);

    notice("Processing barcode and feature filtering options....");
    // load the exclude and include barcode lists
    std::set<std::string> include_id_set;
    std::set<std::string> exclude_id_set;
    if (!include_id_list.empty()) {
        if ( !ignore_ids.empty() ) 
            notice("WARNING: --ignore-ids %s will not take effect if --include-id-list is specified", ignore_ids.c_str());
        tsv_reader id_tr(include_id_list.c_str());
        while (id_tr.read_line() > 0)
            include_id_set.insert(id_tr.str_field_at(0));
    }
    else if (!exclude_id_list.empty()) {
        if ( !ignore_ids.empty() ) 
            notice("WARNING: --ignore-ids %s will not take effect if --exclude-id-list is specified", ignore_ids.c_str());
        tsv_reader id_tr(exclude_id_list.c_str());
        while (id_tr.read_line() > 0)
            exclude_id_set.insert(id_tr.str_field_at(0));
    }
    else {
        // parse ignore_ids
        std::vector<std::string> ignore_ids_vec;
        split(ignore_ids_vec, ",", ignore_ids);
        for(int32_t i=0; i < (int32_t)ignore_ids_vec.size(); ++i)
            exclude_id_set.insert(ignore_ids_vec[i]);
    }

    // load the exclude and include gene lists
    std::set<std::string> include_ftr_set;
    std::set<std::string> exclude_ftr_set;
    if (!include_ftr_list.empty()) {
        tsv_reader ftr_tr(include_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
        include_ftr_set.insert(ftr_tr.str_field_at(0));
    }
    else if (!exclude_ftr_list.empty()) {
        tsv_reader ftr_tr(exclude_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
        exclude_ftr_set.insert(ftr_tr.str_field_at(0));
    }

    std::map<std::string, int32_t> ftr2idx;
    std::vector<std::string> ftr_ids;
    std::vector<uint64_t> ftr_cnts;
    std::map<std::string, int32_t> cell2idx;
    std::vector<std::string> cell_ids;
    std::vector<double> cell_x;
    std::vector<double> cell_y;
    std::vector<int32_t> cell_cnts;
    std::map<int32_t, std::map<int32_t, int32_t> > cell2ftr2cnts;

    notice("Reading the input file... %s", pixelf.c_str());
    tsv_reader pix_tr(pixelf.c_str());
    // detect header columns
    int32_t icol_id = -1;
    int32_t icol_ftr = -1;
    int32_t icol_cnt = -1;
    int32_t icol_x = -1;
    int32_t icol_y = -1;
    if ( no_header ) {
        icol_x = idx_col_x - 1;
        icol_y = idx_col_y - 1;
        icol_ftr = idx_col_ftr - 1;
        icol_cnt = idx_col_cnt - 1;
        icol_id = idx_col_id - 1;
    }
    else {
        if ( pix_tr.read_line() < 0 ) {
            error("Cannot read the input file %s", pixelf.c_str());
        }
        for(int32_t i=0; i < pix_tr.nfields; ++i) {
            if ( strcmp(pix_tr.str_field_at(i), in_col_id.c_str()) == 0 )
                icol_id = i;
            else if ( strcmp(pix_tr.str_field_at(i), in_col_ftr.c_str()) == 0 )
                icol_ftr = i;
            else if ( strcmp(pix_tr.str_field_at(i), in_col_cnt.c_str()) == 0 )
                icol_cnt = i;
            else if ( strcmp(pix_tr.str_field_at(i), in_col_x.c_str()) == 0 )
                icol_x = i;
            else if ( strcmp(pix_tr.str_field_at(i), in_col_y.c_str()) == 0 )
                icol_y = i;
        }
        if ( icol_id < 0 )
            error("Cannot find the column %s in the input file", in_col_id.c_str());
        if ( icol_ftr < 0 )
            error("Cannot find the column %s in the input file", in_col_ftr.c_str());
        if ( icol_cnt < 0 )
            error("Cannot find the column %s in the input file", in_col_cnt.c_str());
    }
    if ( compute_xy ) {
        if ( icol_x < 0 )
            error("Cannot find the column %s in the input file", in_col_x.c_str());
        if ( icol_y < 0 )
            error("Cannot find the column %s in the input file", in_col_y.c_str());
    }

    uint64_t nskip = 0, npass = 0;
    while( pix_tr.read_line() > 0 ) {
        std::string id = pix_tr.str_field_at(icol_id);
        std::string ftr = pix_tr.str_field_at(icol_ftr);
        int32_t cnt = pix_tr.int_field_at(icol_cnt);

        // check if the id is valid
        if ( include_id_set.size() > 0 && include_id_set.find(id) == include_id_set.end() ) {
            nskip++;
            continue;
        }
        if ( exclude_id_set.size() > 0 && exclude_id_set.find(id) != exclude_id_set.end() ) {
            nskip++;
            continue;
        }

        // check if the feature is valid
        if ( include_ftr_set.size() > 0 && include_ftr_set.find(ftr) == include_ftr_set.end() ) {
            nskip++;
            continue;
        }
        if ( exclude_ftr_set.size() > 0 && exclude_ftr_set.find(ftr) != exclude_ftr_set.end() ) {
            nskip++;
            continue;
        }

        npass++;

        std::map<std::string, int32_t>::iterator it_cell = cell2idx.find(id);
        if ( it_cell == cell2idx.end() ) {
            cell2idx[id] = cell_ids.size();
            cell_ids.push_back(id);
            cell_x.push_back(0);
            cell_y.push_back(0);
            cell_cnts.push_back(0);
            it_cell = cell2idx.find(id);
        }
        std::map<std::string, int32_t>::iterator it_ftr = ftr2idx.find(ftr);
        if ( it_ftr == ftr2idx.end() ) {
            ftr2idx[ftr] = ftr_ids.size();
            ftr_ids.push_back(ftr);
            ftr_cnts.push_back(0);
            it_ftr = ftr2idx.find(ftr);
        }
        ftr_cnts[it_ftr->second] += cnt;
        cell2ftr2cnts[it_cell->second][it_ftr->second] += cnt;
        if ( compute_xy ) {
            cell_x[it_cell->second] += (cnt * pix_tr.double_field_at(icol_x));
            cell_y[it_cell->second] += (cnt * pix_tr.double_field_at(icol_y));
            cell_cnts[it_cell->second] += cnt;
        }

        if ( npass % 1000000 == 0 )
            notice("Number of passed lines: %" PRIu64 ", skipped lines: %" PRIu64, npass, nskip);
    }
    notice("Number of skipped lines: %" PRIu64, nskip);
    notice("Number of passed lines: %" PRIu64, npass);

    // sort features based on the names, after thresholding by the min_feature_count
    std::map<std::string, int32_t> ftr2sortidx;
    std::vector<int32_t> sortidx;
    std::vector<int32_t> ftr_idx2sortidx(ftr_ids.size(), -1);
    int32_t isort = 0;
    for(int32_t i=0; i < ftr_ids.size(); ++i) {
        if ( ftr_cnts[i] < min_feature_count )
            continue;
        ftr2sortidx[ftr_ids[i]] = isort;
        sortidx.push_back(i);
        ftr_idx2sortidx[i] = isort;
        isort++;
    }

    // write the sptsv file
    notice("Writing sparse TSV output...");
    char buf[255];
    std::vector<uint64_t> sortidx_sums(sortidx.size(), 0);
    htsFile *wh_tsv = hts_open(out_tsv.c_str(), out_tsv.compare(out_tsv.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    htsFile *wh_cell = NULL;
    if ( !skip_cell_tsv ) {
        wh_cell = hts_open(out_cell_meta.c_str(), out_cell_meta.compare(out_cell_meta.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
        hprintf(wh_cell, "cell_id\tX\tY\tgenes\tcounts\n");
    }
    int32_t n_valid_cells = 0;
    for( std::map<int32_t, std::map<int32_t, int32_t> >::iterator it1 = cell2ftr2cnts.begin(); 
         it1 != cell2ftr2cnts.end(); ++it1 ) {
        // calculate cell count
        int32_t cell_cnt = 0;
        std::map<int32_t, int32_t> sortidx2cnt;
        for( std::map<int32_t, int32_t>::iterator it2 = it1->second.begin(); 
             it2 != it1->second.end(); ++it2 ) {
            int32_t isort = ftr_idx2sortidx[it2->first];
            if ( isort < 0 )
                continue;
            sortidx2cnt[isort] += it2->second;
            cell_cnt += it2->second;
        }
        if ( cell_cnt < min_cell_count )
            continue;
        // print the cell
        snprintf(buf, 255, "%08x", (uint32_t)rand() ); // random key
        hprintf(wh_tsv, "%s\t%s", buf, cell_ids[it1->first].c_str());
        if ( add_xy ) {
            double x = cell_x[it1->first] / cell_cnts[it1->first];
            double y = cell_y[it1->first] / cell_cnts[it1->first];
            hprintf(wh_tsv, "\t%.3f\t%.3f", x, y); 
        }
        hprintf(wh_tsv, "\t%zu\t%d", sortidx2cnt.size(), cell_cnt);
        if ( wh_cell != NULL ) {
            double x = cell_x[it1->first] / cell_cnts[it1->first];
            double y = cell_y[it1->first] / cell_cnts[it1->first];
            hprintf(wh_cell, "%s\t%.3f\t%.3f\t%zu\t%d\n", cell_ids[it1->first].c_str(), x, y, sortidx2cnt.size(), cell_cnt);
        }
        for(std::map<int32_t,int32_t>::iterator it2 = sortidx2cnt.begin(); it2 != sortidx2cnt.end(); ++it2) {
            int32_t isort = it2->first;
            int32_t cnt = it2->second;
            sortidx_sums[isort] += cnt;
            hprintf(wh_tsv, "\t%d %d", isort, cnt);
        }
        hprintf(wh_tsv, "\n");
        n_valid_cells++;
        if ( n_valid_cells % 10000 == 0 )
            notice("Writing %d valid cells", n_valid_cells);
    }
    hts_close(wh_tsv);
    notice("Finished writing %d valid cells", n_valid_cells);

    // write the feature counts
    notice("Writing feature counts...");
    htsFile *wh_counts =
        hts_open(out_feature_counts.c_str(),
                out_feature_counts.compare(out_feature_counts.size() - 3, 3,
                                            ".gz", 3) == 0
                    ? "wz"
                    : "w");
    for (int32_t i = 0; i < (int32_t)sortidx_sums.size(); ++i) {
        hprintf(wh_counts, "%s\t%zu\n", ftr_ids[sortidx[i]].c_str(), sortidx_sums[i]);
    }
    hts_close(wh_counts);

    notice("Writing metadata file... %s", out_json.c_str());
    // write the metadata file
    htsFile *wh_json = hts_open(
        out_json.c_str(),
        out_json.compare(out_json.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    nlohmann::json j_meta;
    std::map<std::string, int32_t> j_dict;
    for (int32_t i = 0; i < (int32_t)sortidx.size(); ++i) {
        j_dict[ftr_ids[sortidx[i]]] = i; // 0-based index for the output feature
    }
    std::vector<std::string> j_header_info;
    j_header_info.push_back(out_col_random_key);
    j_header_info.push_back(out_col_id);
    if ( add_xy ) {
        j_header_info.push_back(out_col_x);
        j_header_info.push_back(out_col_y);
    }
    j_meta[keyname_dictionary] = j_dict;
    j_meta[keyname_header_info] = j_header_info;
    j_meta[keyname_n_features] = (int32_t)sortidx.size(); // number of features
    j_meta[keyname_n_modalities] = 1;                   // number of modalities
    j_meta[keyname_offset_data] = add_xy ? 4 : 2;   // offset for the data columns
    j_meta[keyname_n_units] = n_valid_cells; // number of units
    hprintf(wh_json, "%s\n", j_meta.dump(4).c_str());
    hts_close(wh_json);

    notice("Analysis finished");

    return 0;
}
