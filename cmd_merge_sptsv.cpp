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
// merge-sptsv : Merge multiple sparse TSV files into a single sparse TSV file
/////////////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdMergeSpTSV(int32_t argc, char **argv)
{
    std::string in_list; // input list file containing input sparse TSV files
    std::string out_prefix; // output prefix for merged sparse TSV files
    std::string colname_random_key("random_key");
    int32_t min_feature_count = 0;  // minimum feature count to include in the output
    std::string keyname_dictionary("dictionary");
    std::string keyname_header_info("header_info");
    std::string keyname_n_features("n_features");
    std::string keyname_n_units("n_units");
    std::string keyname_n_modalities("n_modalities");
    std::string keyname_offset_data("offset_data");
    std::string colname_sample_id("sample_id");
    int32_t min_feature_count_per_sample = 0; // minimum feature count per sample to include in the output
    bool exclude_random_key = false;
    bool skip_sample_id_column = false;
    bool combine_headers_to_id = false;
    std::string delimiter(":");
    std::string colname_combined_id("cell_id");
    std::string suffix_feature_counts(".feature.counts.tsv");
    std::string suffix_sptsv(".tsv");
    std::string suffix_json(".json");

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("list", &in_list, "Input list file containing input sparse TSV files")
    LONG_STRING_PARAM("out", &out_prefix, "Output model file to store (gene x factor) matrix")

    LONG_PARAM_GROUP("Header Merging Options", NULL)
    LONG_PARAM("combine-headers-to-id", &combine_headers_to_id, "Combine multiple header fields to create a combined sample ID, when the headers are inconsistent across input files")
    LONG_STRING_PARAM("delim", &delimiter, "Delimiter to use as barcode names when multiple headers fields are combined")
    LONG_STRING_PARAM("colname-combined-id", &colname_combined_id, "Column name for the combined cell ID in the output")
    LONG_STRING_PARAM("colname-sample-id", &colname_sample_id, "Column name for the sample ID in the output")
    LONG_PARAM("skip-sample-id-column", &skip_sample_id_column, "Skip the sample ID column in the output")
    LONG_INT_PARAM("min-feature-count-per-sample", &min_feature_count_per_sample, "Minimum feature count per sample to include in the output.")

    LONG_PARAM_GROUP("Auxiliary Input/Output Options", NULL)
    LONG_STRING_PARAM("colname-random-key", &colname_random_key, "Column name for the random key in the output")
    LONG_STRING_PARAM("keyname-dictionary", &keyname_dictionary, "Key name for the dictionary in the metadata file")
    LONG_STRING_PARAM("keyname-header-info", &keyname_header_info, "Key name for the header information in the metadata file")
    LONG_STRING_PARAM("keyname-n-features", &keyname_n_features, "Key name for the number of features in the metadata file")
    LONG_STRING_PARAM("keyname-n-units", &keyname_n_units, "Key name for the number of units in the metadata file")
    LONG_STRING_PARAM("keyname-offset-data", &keyname_offset_data, "Key name for the offset data in the metadata file")
    LONG_PARAM("exclude-random-key", &exclude_random_key, "Exclude the random key in the output")
    LONG_STRING_PARAM("suffix-feature-counts", &suffix_feature_counts, "Suffix for the per-sample feature count files")
    LONG_STRING_PARAM("suffix-sptsv", &suffix_sptsv, "Suffix for the per-sample sparse TSV files")
    LONG_STRING_PARAM("suffix-json", &suffix_json, "Suffix for the per-sample JSON metadata files")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    // read the list file
    std::vector<std::string> sample_ids;
    std::vector<std::string> feature_count_files;
    std::vector<std::string> sptsv_files;
    std::vector<std::string> json_files;

    tsv_reader tr_list(in_list.c_str());
    while( tr_list.read_line() ) {
        if ( tr_list.nfields == 2 ) {
            const char* sample_id = tr_list.str_field_at(0);
            std::string prefix(tr_list.str_field_at(1));
            sample_ids.push_back(sample_id);
            feature_count_files.push_back(prefix + suffix_feature_counts);
            sptsv_files.push_back(prefix + suffix_sptsv);
            json_files.push_back(prefix + suffix_json);
        }
        else if ( tr_list.nfields == 4 ) {
            sample_ids.push_back(tr_list.str_field_at(0));
            feature_count_files.push_back(tr_list.str_field_at(1));
            sptsv_files.push_back(tr_list.str_field_at(2));
            json_files.push_back(tr_list.str_field_at(3));
        }
        else {
            error("Input list file %s must have 2 or 4 columns: sample_id <prefix> OR sample_id <feature_count_file> <sptsv_file> <json_file>", in_list.c_str());
        }    
    }

    // read the feature count files and determine the features to include
    std::set<std::string> features_to_include;
    std::map<std::string, std::vector<int32_t> > feature_sample_counts; // feature to per-sample counts
    int32_t n_samples = (int32_t)sample_ids.size();
    if ( min_feature_count_per_sample > 0 ) {
        notice("Determining features to include based on minimum feature count per sample %d", min_feature_count_per_sample);
        std::map<std::string, int32_t> ftr2pass; // feature to number of samples passing the threshold
        for(int32_t i=0; i < n_samples; ++i) {
            tsv_reader tr_ftr(feature_count_files[i].c_str());
            while( tr_ftr.read_line() ) {
                std::string ftr_name(tr_ftr.str_field_at(0));
                int32_t cnt = tr_ftr.int_field_at(1);
                if ( cnt >= min_feature_count_per_sample ) {
                    ftr2pass[ftr_name] += 1;
                }
                if ( feature_sample_counts.find(ftr_name) == feature_sample_counts.end() ) {
                    feature_sample_counts[ftr_name].resize(n_samples, 0);
                }
                feature_sample_counts[ftr_name][i] = cnt;
            }
        }
        for( std::map<std::string, int32_t>::const_iterator it = ftr2pass.begin(); it != ftr2pass.end(); ++it ) {
            if ( it->second == n_samples ) {
                features_to_include.insert(it->first);
            }
        }
        notice("Total %d features to be included in the output", (int32_t)features_to_include.size());
    }

    // read JSON files from 
    std::vector<nlohmann::json> jsons(n_samples);
    for(int32_t i=0; i < n_samples; ++i) {
        std::ifstream json_file(json_files[i]);
        if (!json_file.is_open()) {
            error("Cannot open metadata file %s", json_files[i].c_str());
        }
        json_file >> jsons[i];
    }

    std::vector<int32_t> icols_random_keys(n_samples, -1);
    std::vector< std::vector<std::string> > header_keys(n_samples);
    std::vector<int32_t> in_offset_data(n_samples, -1);
    std::vector<int32_t> in_n_features(n_samples, -1);
    std::vector<int32_t> in_n_units(n_samples, -1);
    // merge the header information
    for(int32_t i=0; i < n_samples; ++i) {
        try {
            std::vector<std::string> hdrs = jsons[i][keyname_header_info].get<std::vector<std::string> >();
            header_keys[i] = hdrs;
            for(int32_t j=0; j < (int32_t)hdrs.size(); ++j) {
                if ( hdrs[j].compare(colname_random_key) == 0 ) {
                    icols_random_keys[i] = j;
                }
            }
            in_offset_data[i] = jsons[i][keyname_offset_data].get<int32_t>();
            in_n_features[i] = jsons[i][keyname_n_features].get<int32_t>();
            in_n_units[i] = jsons[i][keyname_n_units].get<int32_t>();
        }
        catch (nlohmann::json::exception& e) {
            error("Cannot read header information from the metadata file %s: %s", json_files[i].c_str(), e.what());
        }
        // check whether the header keys are consistent
    }

    for(int32_t i=1; i < n_samples; ++i) {
        bool consistent = true;
        if ( header_keys[i].size() != header_keys[0].size() ) {
            consistent = false;
        }
        else {
            for(int32_t j=0; j < (int32_t)header_keys[0].size(); ++j) {
                if ( header_keys[i][j].compare( header_keys[0][j] ) != 0 ) {
                    consistent = false;
                    break;
                }
            }
        }
        if ( !consistent ) {
            if ( !combine_headers_to_id ) {
                error("Header information is inconsistent across input files. Use --combine-headers-to-id to combine multiple header fields to create a unique cell ID.");
            }
            else {
                notice("Header information is inconsistent across input files, but they will be combined to create unique IDs.");
                break;
            }
        }
    }
    
    // std::vector<std::string> merged_header_keys;
    // if ( combine_headers_to_id ) {
    //     for(int32_t i=0; i < n_samples; ++i) {
    //         std::vector<std::string> hdrs;
    //         if ( !skip_sample_id_column ) {
    //             hdrs.push_back(sample_ids[i]);
    //         }
    //         for(int32_t j=0; j < (int32_t)header_keys[i].size(); ++j) {
    //             if ( j != icols_random_keys[i] ) { // skip the random key column
    //                 hdrs.push_back( header_keys[i][j] );
    //             }
    //         }
    //         std::string combined_id(hdrs[0]);
    //         for(int32_t j=1; j < (int32_t)hdrs.size(); ++j) {
    //             combined_id += (delimiter + hdrs[j]);
    //         }
    //         merged_header_keys.push_back(combined_id);
    //     }
    // }

    std::vector<std::string> out_header_info;
    std::vector< std::vector<int32_t> > out_header_key_indices(n_samples);
    std::vector< std::vector<int32_t> > out_header_combined_indices(n_samples);
    if ( !exclude_random_key ) {
        out_header_info.push_back(colname_random_key);
        for(int32_t i=0; i < n_samples; ++i) {
            out_header_key_indices[i].push_back( icols_random_keys[i] );
        }
    }
    if ( !skip_sample_id_column ) {
        out_header_info.push_back(colname_sample_id);
        for(int32_t i=0; i < n_samples; ++i) {
            out_header_key_indices[i].push_back( -1 ); // sample ID column
        }
    }
    if ( combine_headers_to_id ) {
        out_header_info.push_back(colname_combined_id);
        for(int32_t i=0; i < n_samples; ++i) {
            out_header_key_indices[i].push_back( -2 ); // combined ID column
            if ( skip_sample_id_column) {
                out_header_combined_indices[i].push_back( -1 ); // sample ID column
            }
            for(int32_t j=0; j < (int32_t)header_keys[i].size(); ++j) {
                if ( j != icols_random_keys[i] ) { // skip the random key column
                    out_header_combined_indices[i].push_back( j );
                }
            }
        }
    }
    else {
        for(int32_t i=0; i < n_samples; ++i) {
            for(int32_t j=0; j < (int32_t)header_keys[i].size(); ++j) {
                if ( j != icols_random_keys[i] ) { // skip the random key column
                    if ( i == 0 ) {
                        out_header_info.push_back( header_keys[i][j] );
                    }
                    out_header_key_indices[i].push_back( j );
                }
            }
        }
    }
    int32_t out_offset_data = out_header_info.size();

    // create merged feature dictionary
    std::vector< std::map<int32_t, int32_t> > in2out_feature_indices(n_samples);
    std::vector< std::map<std::string, int32_t> > in_feature_name2idx(n_samples);
    // construct union of all features
    std::map<std::string, int32_t> out_feature_name2idx;
    for(int32_t i=0; i < n_samples; ++i) {
        try {
            in_feature_name2idx[i] = jsons[i][keyname_dictionary].get<std::map<std::string, int32_t> >();
            for (const auto& kv : in_feature_name2idx[i]) {
                // check if the feature is to be included
                if ( ( min_feature_count_per_sample == 0 ) ||
                     ( features_to_include.find(kv.first) != features_to_include.end() ) ) {
                    std::map<std::string, int32_t>::iterator it = out_feature_name2idx.find(kv.first);
                    if ( it == out_feature_name2idx.end() ) {
                        // add new feature
                        int32_t out_idx = (int32_t)out_feature_name2idx.size();
                        out_feature_name2idx[kv.first] = out_idx;
                        in2out_feature_indices[i][kv.second] = out_idx;
                    }
                    else {
                        // use existing feature
                        in2out_feature_indices[i][kv.second] = it->second;
                    }
                }
            }
        }
        catch (nlohmann::json::exception& e) {
            error("Cannot read dictionary from the metadata file %s: %s", json_files[i].c_str(), e.what());
        }
    }

    // write the merged sparse TSV file
    std::string out_sptsv_file = out_prefix + suffix_sptsv;
    htsFile* wf = hts_open(out_sptsv_file.c_str(), out_sptsv_file.compare(out_sptsv_file.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "wt");
    std::map<int32_t, int32_t> out_feature_counts;
    if ( wf == NULL ) {
        error("Cannot open output sparse TSV file %s for writing", out_sptsv_file.c_str());
    }
    int32_t n_units = 0;
    for(int32_t i=0; i < n_samples; ++i) {
        tsv_reader tr(sptsv_files[i].c_str());
        std::vector<int32_t>& hdr_indices = out_header_key_indices[i];
        while( tr.read_line() ) {
            // print header fields;
            std::string combined_id;
            if ( combine_headers_to_id ) { // construct combined ID
                std::vector<int32_t>& indices = out_header_combined_indices[i];
                for(int32_t j=0; j < (int32_t)indices.size(); ++j) {
                    int32_t hdr_idx = indices[j];
                    if ( j > 0 ) {
                        combined_id += delimiter;
                    }
                    if ( hdr_idx == -1 ) {
                        combined_id += sample_ids[i];
                    }
                    else {
                        combined_id += tr.str_field_at(hdr_idx);
                    }
                }
            }
            for(int32_t j=0; j < (int32_t)hdr_indices.size(); ++j) {
                if ( j > 0 ) {
                    hprintf(wf, "\t");
                }
                int32_t hdr_idx = hdr_indices[j];
                if ( hdr_idx == -1 ) {
                    // sample ID column
                    hprintf(wf, "%s", sample_ids[i].c_str());
                }
                else if ( hdr_idx == -2 ) {
                    // combined ID column
                    hprintf(wf, "%s", combined_id.c_str());
                }
                else {
                    hprintf(wf, "%s", tr.str_field_at(hdr_idx));
                }
            }
            // print the rest of data
            int32_t offset_data = in_offset_data[i];
            std::vector<int32_t> feature_indices;
            std::vector<int32_t> feature_counts;
            int32_t n_genes = tr.int_field_at(offset_data);
            int32_t n_mols = tr.int_field_at(offset_data + 1);
            if ( tr.nfields != 2 * n_genes + offset_data + 2 ) {
                error("The number of fields in the input TSV file (%d) does not match the expected number (%d)", tr.nfields, n_genes + offset_data + 2);
            }
            int32_t sum_cnt = 0;
            for(int32_t j = 0; j < n_genes; ++j) {
                int32_t iftr = tr.int_field_at(2*j + offset_data + 2);
                int32_t cnt = tr.int_field_at(2*j + offset_data + 3);
                std::map<int32_t, int32_t>::iterator it = in2out_feature_indices[i].find(iftr);
                if ( it != in2out_feature_indices[i].end() ) {
                    // feature to be included
                    feature_indices.push_back( it->second );
                    feature_counts.push_back( cnt );
                    sum_cnt += cnt;
                    out_feature_counts[ it->second ] += cnt;
                }
            }
            // write down the contents
            hprintf(wf, "\t%zu\t%d", feature_indices.size(), sum_cnt);
            for(int32_t j=0; j < (int32_t)feature_indices.size(); ++j) {
                hprintf(wf, "\t%d %d", feature_indices[j], feature_counts[j]);
            }
            hprintf(wf, "\n");
            ++n_units;
        }
    }
    hts_close(wf);

    // write down the merged feature count file
    notice("Writing output feature count file...");
    std::string out_feature_counts_file = out_prefix + suffix_feature_counts;
    wf = hts_open(out_feature_counts_file.c_str(), out_feature_counts_file.compare(out_feature_counts_file.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "wt");
    if ( wf == NULL ) {
        error("Cannot open output feature count file %s for writing", out_sptsv_file.c_str());
    }
    for( std::map<std::string, int32_t>::const_iterator it = out_feature_name2idx.begin(); it != out_feature_name2idx.end(); ++it ) {
        int32_t ftr_idx = it->second;
        int32_t cnt = out_feature_counts[ftr_idx];
        hprintf(wf, "%s\t%d\n", it->first.c_str(), cnt);
    }
    hts_close(wf);
    
    // write down the merged JSON metadata file
    notice("Writing output metadata file...");
    std::string out_json_file = out_prefix + suffix_json;
    nlohmann::json out_json;
    out_json[keyname_dictionary] = out_feature_name2idx;
    out_json[keyname_header_info] = out_header_info;
    out_json[keyname_n_features] = (int32_t)out_feature_name2idx.size();
    out_json[keyname_n_units] = n_units;
    out_json[keyname_n_modalities] = 1;                   // number of modalities
    out_json[keyname_offset_data] = out_offset_data;
    htsFile *wh_json = hts_open(
        out_json_file.c_str(),
        out_json_file.compare(out_json_file.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wh_json == NULL ) {
        error("Cannot open output metadata file %s for writing", out_json_file.c_str());
    }
    hprintf(wh_json, "%s\n", out_json.dump(4).c_str());
    hts_close(wh_json);

    notice("Analysis finished");

    return 0;
}
