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


/////////////////////////////////////////////////////////////////////////////////////////////////
// sptsv2model : Create model matrix from Sparse TSV format in FICTURE2 with cluster assignment
/////////////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdSpTSV2Model(int32_t argc, char **argv)
{
    std::string in_tsv;
    std::string in_meta;
    std::string in_features;
    std::string in_clust;
    std::string in_fit;
    std::string out_model;
    std::string colname_random_key("random_key");
    int32_t min_feature_count = 0;  // minimum feature count to include in the output
    std::string keyname_dictionary("dictionary");
    std::string keyname_header_info("header_info");
    std::string keyname_n_features("n_features");
    std::string keyname_n_units("n_units");
    std::string keyname_offset_data("offset_data");
    bool include_random_key = false;
    std::string delim(":");

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input Options", NULL)
    LONG_STRING_PARAM("tsv", &in_tsv, "Input TSV file containing the sparse matrix data")
    LONG_STRING_PARAM("json", &in_meta, "Input JSON metadata file containing the header information")
    LONG_STRING_PARAM("features", &in_features, "Input features file containing the feature names and counts")
    LONG_STRING_PARAM("clust", &in_clust, "Input cluster file containing the cluster assignments for each barcode")
    LONG_STRING_PARAM("fit", &in_fit, "Input fit results file containing the probabilistic cluster assignment for each barcode")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_STRING_PARAM("out", &out_model, "Output model file to store (gene x factor) matrix")

    LONG_PARAM_GROUP("Auxiliary Input/Output Options", NULL)
    LONG_INT_PARAM("min-count", &min_feature_count, "Minimum feature count to include in the output.")
    LONG_STRING_PARAM("colname-random-key", &colname_random_key, "Column name for the random key in the output")
    LONG_PARAM("include-random-key", &include_random_key, "Include the random key in the output")
    LONG_STRING_PARAM("keyname-dictionary", &keyname_dictionary, "Key name for the dictionary in the metadata file")
    LONG_STRING_PARAM("keyname-header-info", &keyname_header_info, "Key name for the header information in the metadata file")
    LONG_STRING_PARAM("keyname-n-features", &keyname_n_features, "Key name for the number of features in the metadata file")
    LONG_STRING_PARAM("keyname-n-units", &keyname_n_units, "Key name for the number of units in the metadata file")
    LONG_STRING_PARAM("keyname-offset-data", &keyname_offset_data, "Key name for the offset data in the metadata file")
    LONG_STRING_PARAM("delim", &delim, "Delimiter to use as barcode names when multiple headers fields are combined")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if ( out_model.empty() ) {
        error("--out must be specified");
    }

    if ( in_tsv.empty() || in_meta.empty() ) {
        error("--tsv, --json, must be specified");
    }

    if ( in_clust.empty() && in_fit.empty() ) {
        error("Either --clust and --fit must be specified");
    }

    if ( !in_clust.empty() && !in_fit.empty() ) {
        notice("Both --clust and --fit are specified. Use only one option");
    }

    if ( in_features.empty() && min_feature_count > 0 ) {
        error("--features must be specified if --min-count is greater than 0");
    }

    // Read the feature count file
    std::map<std::string, int32_t> ftr2cnt;
    if (!in_features.empty()) {
        tsv_reader tr_ftr(in_features.c_str());
        while( tr_ftr.read_line() ) {
            std::string ftr_name = tr_ftr.str_field_at(0);
            if ( min_feature_count > 0 ) {
                if ( tr_ftr.nfields < 2 ) {
                    error("Feature file %s does not contain count, but minimum feature count is set to %d", in_features.c_str(), min_feature_count);
                }
                int32_t cnt = tr_ftr.int_field_at(1);
                if ( cnt >= min_feature_count ) {
                    ftr2cnt[ftr_name] = cnt;
                }
            }
            else {
                int cnt = tr_ftr.nfields < 2 ? 0 : tr_ftr.int_field_at(1);
                ftr2cnt[ftr_name] = cnt; // store the feature count
            }
        }
        tr_ftr.close();
    }

    // Read the metadata file
    notice("Reading metadata file %s", in_meta.c_str());
    nlohmann::json json_data;
    int32_t icol_random_key = -1; // column index for the random key, -1 if not present
    int32_t offset_data = -1;     // default offset for the data in the input TSV file
    int32_t n_features = -1; // number of features in the input TSV file
    int32_t n_units = -1;    // number of units in the input TSV file
    std::vector<std::string> feature_names;
    std::map<std::string, int32_t> ftr2idx; 
    {
        std::ifstream meta_file(in_meta);
        if (!meta_file.is_open()) {
            error("Cannot open metadata file %s", in_meta.c_str());
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
            error("Cannot read header information from the metadata file %s: %s", in_meta.c_str(), e.what());
        }
        offset_data = json_data["offset_data"].get<int32_t>();
        n_features = json_data[keyname_n_features].get<int32_t>();
        n_units = json_data[keyname_n_units].get<int32_t>();
        try {
            std::map<std::string, int32_t> dict = json_data[keyname_dictionary].get<std::map<std::string, int32_t> >();
            feature_names.resize(n_features);
            for (const auto& kv : dict) {
                if ( kv.second < 0 || kv.second >= n_features ) {
                    error("Invalid feature index %d for feature %s in the metadata file %s", kv.second, kv.first.c_str(), in_meta.c_str());
                }
                feature_names[kv.second] = kv.first; // map feature index to feature name
                ftr2idx[kv.first] = kv.second; // store the feature index
            }
            // check if all features are present
            for(int32_t i = 0; i < n_features; ++i) {
                if ( feature_names[i].empty() ) {
                    error("Feature %d is not present in the metadata file %s", i, in_meta.c_str());
                }
            }
        }
        catch (nlohmann::json::exception& e) {
            error("Cannot read dictionary from the metadata file %s: %s", in_meta.c_str(), e.what());
        }

        if ( icol_random_key < 0 && include_random_key ) {
            error("Random key column %s not found in the metadata file %s", colname_random_key.c_str(), in_meta.c_str());      
        }
        if ( offset_data < 0 ) {
            error("Invalid offset_data value %d in the metadata file %s", offset_data, in_meta.c_str());
        }
        if ( n_features <= 0 ) {
            error("Invalid number of features %d in the metadata file %s", n_features, in_meta.c_str());
        }
        if ( n_units <= 0 ) {
            error("Invalid number of units %d in the metadata file %s", n_units, in_meta.c_str());
        }
    }

    // read the cluster assignment information
    std::vector<std::string> clust_ids;
    std::map<std::string, int32_t> clust2idx; // map cluster ID to index
    std::map<std::string, int32_t> unit2idx;
    std::map<std::string, std::vector<double> > unit2probs;
    std::vector<std::string> sorted_clust_ids;
    std::vector<int32_t> srt2unsrt_idx;
    bool mode_clust = false;

    if ( !in_clust.empty() ) {     
        notice("Reading cluster assignment file %s", in_clust.c_str());
        mode_clust = true; // we are in cluster mode
        tsv_reader tsv_tr(in_clust.c_str());
        if ( !tsv_tr.open(in_clust.c_str()) ) {
            error("Cannot open input TSV file %s for reading", in_clust.c_str());
        }   
        int32_t offset_tsv = offset_data + (include_random_key ? 0 : -1);
        while( tsv_tr.read_line() ) {
            std::string barcode;
            for(int32_t i = 0; i < offset_tsv; ++i) {
                if ( barcode.empty() ) {
                    barcode.assign(tsv_tr.str_field_at(i));
                }
                else {
                    barcode += (delim + tsv_tr.str_field_at(i)); // combine multiple fields with the specified delimiter
                }
            }

            const char* clust = tsv_tr.str_field_at(offset_tsv);
            std::map<std::string, int32_t>::iterator it = clust2idx.find(clust);
            int32_t idx = -1;
            if ( it == clust2idx.end() ) {
                clust_ids.push_back(clust);
                idx = (int32_t)clust_ids.size() - 1;
            }
            else {
                idx = it->second;
            }
            clust2idx[clust] = idx; // map cluster ID to index
            unit2idx[barcode] = idx;
        }

        // sort the clusters as numbers, and alphanumerically if not numbers
        sorted_clust_ids = clust_ids;
        std::sort(sorted_clust_ids.begin(), sorted_clust_ids.end(), [](const std::string& a, const std::string& b) {
            if ( a.find_first_not_of("0123456789") == std::string::npos && b.find_first_not_of("0123456789") == std::string::npos ) {
                return std::stoi(a) < std::stoi(b);
            }
            else {
                return a < b;
            }
        });

        std::map<std::string, int32_t> clust2unsrt_idx;
        for(int32_t i = 0; i < (int32_t)clust_ids.size(); ++i) {
            clust2unsrt_idx[clust_ids[i]] = i; // map cluster ID to sorted index
        }
        for(int32_t i = 0; i < (int32_t)sorted_clust_ids.size(); ++i) {
            srt2unsrt_idx.push_back(clust2unsrt_idx[sorted_clust_ids[i]]);
        }
        // for(int32_t i = 0; i < (int32_t)clust_ids.size(); ++i) {
        //     notice("%d\t%s\t%d\t%s", i, clust_ids[i].c_str(), srt2unsrt_idx[i], sorted_clust_ids[i].c_str());
        // }
    }
    else { // in_fit is not empty
        notice("Reading fit results file %s", in_fit.c_str());
        tsv_reader tsv_tr(in_fit.c_str());
        if ( !tsv_tr.open(in_fit.c_str()) ) {
            error("Cannot open input TSV file %s for reading", in_clust.c_str());
        }   
        uint64_t n_lines = 0;
        int32_t offset_tsv = offset_data + (include_random_key ? 0 : -1);
        while( tsv_tr.read_line() ) {
            if ( n_lines == 0 ) {
                // read the header line to get the cluster IDs
                for(int32_t i = offset_tsv; i < tsv_tr.nfields; ++i) {
                    const char* id = tsv_tr.str_field_at(i);
                    clust_ids.push_back(id);
                    clust2idx[id] = (int32_t)clust_ids.size() - 1; // map cluster ID to index
                }    
            }
            else {
                std::string barcode;
                for(int32_t i = 0; i < offset_tsv; ++i) {
                    if ( barcode.empty() ) {
                        barcode.assign(tsv_tr.str_field_at(i));
                    }
                    else {
                        barcode += (delim + tsv_tr.str_field_at(i)); // combine multiple fields with the specified delimiter
                    }
                }
                std::vector<double>& probs = unit2probs[barcode];
                probs.resize(clust_ids.size(), 0.0);
                if ( tsv_tr.nfields != offset_tsv + clust_ids.size() ) {
                    error("The number of fields in the input TSV file (%d) does not match the expected number (%d)", tsv_tr.nfields, offset_tsv + clust_ids.size());
                }
                double sum_probs = 0.0;
                for(int32_t i = 0; i < (int32_t)clust_ids.size(); ++i) {
                    probs[i] = tsv_tr.double_field_at(offset_tsv + i); // store the probabilities for each cluster
                    sum_probs += probs[i];
                }
                if ( fabs(sum_probs - 1.0) > 1e-5 ) { // normalization needed
                    for(int32_t i = 0; i < (int32_t)clust_ids.size(); ++i) {
                        probs[i] = ( probs[i] + 1e-20 ) / ( sum_probs + 1e-20 * clust_ids.size() ); // normalize the probabilities to sum to 1
                    }
                }
            }
            ++n_lines;
        }

        // do not sort the cluster in probability mode
        sorted_clust_ids = clust_ids; // use the original order of cluster IDs
        for(int32_t i = 0; i < (int32_t)clust_ids.size(); ++i) {
            srt2unsrt_idx.push_back(i);
        }
    }

    // read the input TSV file and create the model matrix
    std::map<int32_t, std::vector<double> > iftr2vec; // feature to vector of counts
    std::map<int32_t, std::vector<double> >::iterator iftr2vec_it;
    std::map<std::string, int32_t>::iterator unit2idx_it;
    std::map<std::string, std::vector<double> >::iterator unit2probs_it;
    int32_t n_skip_units = 0;


    tsv_reader tsv_tr(in_tsv.c_str());
    if ( !tsv_tr.open(in_tsv.c_str()) ) {
        error("Cannot open input TSV file %s for reading", in_tsv.c_str());
    }   
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
        int32_t n_genes = tsv_tr.int_field_at(offset_data);
        int32_t n_mols = tsv_tr.int_field_at(offset_data + 1);
        if ( tsv_tr.nfields != 2 * n_genes + offset_data + 2 ) {
            error("The number of fields in the input TSV file (%d) does not match the expected number (%d)", tsv_tr.nfields, n_genes + offset_data + 2);
        }

        // identify the current cluster index
        int32_t clust_idx = -1;
        std::vector<double> clust_probs;
        if ( mode_clust ) {
            unit2idx_it = unit2idx.find(barcode);
            if ( unit2idx_it == unit2idx.end() ) {
                ++n_skip_units;
                if ( n_skip_units < 10 ) {
                    notice("Barcode %s not found in the cluster assignment file %s... Skipping", barcode.c_str(), in_clust.c_str());
                }
                else if ( n_skip_units == 10 ) {
                    notice("Skipping %d barcodes not found in the cluster assignment file %s", n_skip_units, in_clust.c_str());
                }
                continue;
            }
            clust_idx = unit2idx_it->second;
        }
        else {
            unit2probs_it = unit2probs.find(barcode);
            if ( unit2probs_it == unit2probs.end() ) {
                ++n_skip_units;
                if ( n_skip_units < 10 ) {
                    notice("Barcode %s not found in the fit results file %s... Skipping", barcode.c_str(), in_fit.c_str());
                }
                else if ( n_skip_units == 10 ) {
                    notice("Skipping %d barcodes not found in the fit results file %s", n_skip_units, in_fit.c_str());
                }
                continue;
            }
            clust_probs = unit2probs_it->second; 
        }

        int32_t sum_cnt = 0;
        for(int32_t i = 0; i < n_genes; ++i) {
            const char* s = tsv_tr.str_field_at(2*i + offset_data + 2);
            const char* t = tsv_tr.str_field_at(2*i + offset_data + 3);
            if ( t == NULL ) {
                error("Invalid data format in the input TSV file %s at line %d: expected space-separated values", in_tsv.c_str(), tsv_tr.nlines);
            }
            int32_t iftr = atoi(s);
            int32_t cnt = atoi(t);

            iftr2vec_it = iftr2vec.find(iftr);
            if ( iftr2vec_it == iftr2vec.end() ) { // feature not found, create a new vector
                iftr2vec[iftr].resize(clust_ids.size(), 0.0); // initialize the vector for this feature
                iftr2vec_it = iftr2vec.find(iftr);
            }

            if ( mode_clust ) {
                if ( clust_idx >= 0 ) {
                    iftr2vec_it->second[clust_idx] += cnt; // add the count to the cluster index
                }
            }
            else {
                if ( !clust_probs.empty() ) {
                    for(int32_t iclust = 0; iclust < (int32_t)clust_probs.size(); ++iclust) {
                        iftr2vec_it->second[iclust] += (cnt * clust_probs[iclust]); // add the weighted count to the cluster index
                    }
                }
            }
        }
    }
    tsv_tr.close();

    // write the output model file
    notice("Writing output model file %s", out_model.c_str());
    htsFile* wf = hts_open(out_model.c_str(), out_model.compare(out_model.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output model file %s for writing", out_model.c_str());
    }
    hprintf(wf, "Feature");
    for(int32_t i = 0; i < (int32_t)sorted_clust_ids.size(); ++i) {
        hprintf(wf, "\t%s", sorted_clust_ids[i].c_str());
    }
    hprintf(wf, "\n");
    for(iftr2vec_it = iftr2vec.begin(); iftr2vec_it != iftr2vec.end(); ++iftr2vec_it) {
        int32_t iftr = iftr2vec_it->first;
        if ( !ftr2cnt.empty() && ftr2cnt.find(feature_names[iftr]) == ftr2cnt.end() ) {  // skip features not included
            continue; // skip features not in the feature count file
        }
        const std::vector<double>& vec = iftr2vec_it->second;
        hprintf(wf, "%s", feature_names[iftr].c_str());
        for(int32_t i = 0; i < (int32_t)sorted_clust_ids.size(); ++i) {
            hprintf(wf, "\t%.3f", vec[srt2unsrt_idx[i]]);
        }
        hprintf(wf, "\n");
    }
    hts_close(wf);

    notice("Analysis finished");

    return 0;
}
