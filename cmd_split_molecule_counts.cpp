#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

/////////////////////////////////////////////////////////////////////////////////////////
// split-molecule-counts : Split molecule counts into pixel-level bins
/////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdSplitMoleculeCounts(int32_t argc, char **argv)
{
    std::string in_mol_tsv;   // TSV file containing individual molecules
    std::string in_ftr_tsv;  // TSV file containing gene-level total counts (not used in the current version but can be used for filtering lowly expressed genes in the future)
    std::string out_prefix;   // Output Prefix
    std::string in_mol_tsv_delim = "\t"; // Delimiter for the input molecule TSV file
    std::string in_ftr_tsv_delim = "\t"; // Delimiter for the input feature TSV file
    std::string out_mol_tsv_delim = "\t"; // Delimiter for the output molecule TSV file
    std::string out_ftr_tsv_delim = "\t"; // Delimiter for the output feature TSV file
    std::string out_mol_suffix = "molecules.tsv.gz"; // Suffix for the output molecule TSV file
    std::string out_ftr_suffix = "features.tsv.gz"; // Suffix for the output feature TSV file
    std::string out_index_suffix = "_index.tsv"; 
    std::string out_json_suffix = "_bin_counts.json"; // Suffix for the output JSON file containing bin counts and other metadata
    std::string colname_x = "X"; // Column name for X coordinate
    std::string colname_y = "Y"; // Column name for Y coordinate
    std::string colname_feature = "gene"; // Column name for gene name
    std::string colname_count = "count"; // Column name for gene count
    std::string strip_comment_char = "#"; // Character to strip from the beginning of lines in the input files (if any)
    std::vector<std::string> col_renames; // Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...
    int32_t bin_count = 50; // When --equal-bins is used, determine the number of bins to split the data into the same bin
    bool compact_bin = false; // Whether to compact the bins to store minimal information.
    bool skip_original = false; // Whether to skip writing the original file


    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("mol-tsv", &in_mol_tsv, "TSV file containing individual molecules")
    LONG_STRING_PARAM("feature-tsv", &in_ftr_tsv, "TSV file containing gene-level total counts (not used in the current version but can be used for filtering lowly expressed genes in the future)")
    LONG_STRING_PARAM("out-prefix", &out_prefix, "Output Prefix for the joined TSV files")
  
    LONG_PARAM_GROUP("Key Parameters", NULL)
    LONG_INT_PARAM("bin-count", &bin_count, "When --equal-bins is used, determine the number of bins to split the data into the same bin")
    LONG_PARAM("skip-original", &skip_original, "Whether to skip writing the original file")
    LONG_PARAM("compact-bin", &compact_bin, "Compact the bins to store minimal information")

    LONG_PARAM_GROUP("Expected columns in input and output", NULL)
    LONG_STRING_PARAM("colname-feature", &colname_feature, "Column name for gene name")
    LONG_STRING_PARAM("colname-count", &colname_count, "Column name for gene count")
    LONG_STRING_PARAM("colname-x", &colname_x, "Column name for X coordinate")
    LONG_STRING_PARAM("colname-y", &colname_y, "Column name for Y coordinate")
    LONG_MULTI_STRING_PARAM("col-rename", &col_renames, "Columns to rename in the output file. Format: old_name1:new_name1 old_name2:new_name2 ...")

    LONG_PARAM_GROUP("Auxilary Input/Output Parameters", NULL)
    LONG_STRING_PARAM("in-mol-tsv-delim", &in_mol_tsv_delim, "Delimiter for the input molecule TSV file")
    LONG_STRING_PARAM("in-feature-tsv-delim", &in_ftr_tsv_delim, "Delimiter for the input feature TSV file")
    LONG_STRING_PARAM("out-mol-tsv-delim", &out_mol_tsv_delim, "Delimiter for the output molecule TSV file")
    LONG_STRING_PARAM("out-feature-tsv-delim", &out_ftr_tsv_delim, "Delimiter for the output feature TSV file")
    LONG_STRING_PARAM("out-mol-suffix", &out_mol_suffix, "Suffix for the output molecule TSV file")
    LONG_STRING_PARAM("out-json-suffix", &out_json_suffix, "Suffix for the output JSON file containing bin counts and other metadata")
    LONG_STRING_PARAM("out-feature-suffix", &out_ftr_suffix, "Suffix for the output feature TSV file")
    LONG_STRING_PARAM("strip-comment-char", &strip_comment_char, "Character to strip from the beginning of lines in the input files (if any)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    // read the feature count TSV file and assign each gene into ordered bins
    notice("Reading feature count TSV file %s", in_ftr_tsv.c_str());
    tsv_reader tr_ftr(in_ftr_tsv.c_str());
    tr_ftr.delimiter = in_ftr_tsv_delim[0];
    std::map<std::string, int32_t> ftr2cnts;
    uint64_t sum_cnt = 0;
    while (tr_ftr.read_line() > 0) {
        if (tr_ftr.str_field_at(0)[0] == strip_comment_char[0]) {
            continue; // skip comment lines
        }
        int32_t cnt = tr_ftr.int_field_at(1);
        if ( cnt > 0 ) {
            ftr2cnts[tr_ftr.str_field_at(0)] = cnt;            
            sum_cnt += cnt;
        }
    }

    std::map<std::string, std::string> col_rename_map;
    for(int32_t i=0; i < col_renames.size(); ++i) {
        std::string rename_str = col_renames[i];
        size_t colon_pos = rename_str.find(':');
        if ( colon_pos == std::string::npos ) {
            error("Invalid column rename format: %s. Expected format: old_name:new_name", rename_str.c_str());
        }
        std::string old_name = rename_str.substr(0, colon_pos);
        std::string new_name = rename_str.substr(colon_pos + 1);
        col_rename_map[old_name] = new_name;
    }

    // sort the genes by their counts and assign bins
    std::vector<std::pair<std::string, int32_t>> ftr_cnt_vec(ftr2cnts.begin(), ftr2cnts.end());
    std::sort(ftr_cnt_vec.begin(), ftr_cnt_vec.end(), [](const std::pair<std::string, int32_t>& a, const std::pair<std::string, int32_t>& b) {
        return a.second > b.second; // sort in descending order of counts
    });

    // pick genes for each bin
    std::map<std::string, int32_t> ftr2bin;
    uint64_t cumsum_cnt = 0;
    int32_t j = 0;
    std::vector<uint64_t> bin_mol_cnts(bin_count, 0);
    std::vector<int32_t> bin_ftr_cnts(bin_count, 0);
    int32_t actual_bin_count = 0;
    for(int32_t i=0; i < bin_count; ++i) {
        // determine threshold - average number of molecules among the remaining bins
        int32_t bin_thres = i == bin_count - 1 ? INT_MAX : (sum_cnt - cumsum_cnt) / (bin_count - i);
        //notice("sum_cnt = %llu, cumsum_cnt = %llu, remaining_cnt = %llu, bin_thres = %d", sum_cnt, cumsum_cnt, sum_cnt - cumsum_cnt, bin_thres);
        if ( bin_thres == 0 ) { // include at least one gene in each bin if there are still genes left
            bin_thres = 1;
        }
        uint64_t current_bin_cnt = 0;
        int32_t current_bin_ngenes = 0;
        while (j < ftr_cnt_vec.size() && current_bin_cnt < bin_thres) {
            ftr2bin[ftr_cnt_vec[j].first] = i;
            current_bin_cnt += ftr_cnt_vec[j].second;
            ++j;
            ++current_bin_ngenes;
        }
        if ( current_bin_ngenes > 0 ) {
            cumsum_cnt += current_bin_cnt;
            bin_mol_cnts[i] = current_bin_cnt;
            bin_ftr_cnts[i] = current_bin_ngenes;
            notice("Bin %d: %d / %d genes, %llu molecules, (%.3f%% cumulative, threshold = %d)", i, current_bin_ngenes, ftr_cnt_vec.size(), current_bin_cnt, (double)cumsum_cnt / sum_cnt * 100, bin_thres);
            ++actual_bin_count;
        }
        else { // no more genes to assign, stop here
            break;
        }
    }
    notice("Total %d bins assigned, %d genes, %llu molecules", actual_bin_count, (int32_t)ftr2bin.size(), sum_cnt);

    // create output files all at once
    bool out_mol_gz = out_mol_suffix.size() > 3 && out_mol_suffix.substr(out_mol_suffix.size() - 3) == ".gz";
    bool out_ftr_gz = out_ftr_suffix.size() > 3 && out_ftr_suffix.substr(out_ftr_suffix.size() - 3) == ".gz";
    bool out_index_gz = out_index_suffix.size() > 3 && out_index_suffix.substr(out_index_suffix.size() - 3) == ".gz";

    htsFile* wf_mol_all = skip_original ? NULL : hts_open((out_prefix + "_all_" + out_mol_suffix).c_str(), out_mol_gz ?  "wz" : "w");
    std::vector<htsFile*> wf_mol_bins(actual_bin_count, NULL);
    for(int32_t i=0; i < actual_bin_count; ++i) {
        wf_mol_bins[i] = hts_open((out_prefix + "_bin" + std::to_string(i+1) + "_" + out_mol_suffix).c_str(), out_mol_gz ?  "wz" : "w");
    }

    // read the molecule TSV file and write to corresponding bin files
    tsv_reader tr_mol(in_mol_tsv.c_str());
    tr_mol.delimiter = in_mol_tsv_delim[0];
    std::string out_line;
    std::string out_compact_line;
    int32_t col_idx_feature = -1;
    int32_t col_idx_x = -1;
    int32_t col_idx_y = -1;
    int32_t col_idx_count = -1;
    if ( tr_mol.read_line() > 0 ) {
        const char* first_col = tr_mol.str_field_at(0);
        while( strip_comment_char.size() > 0 && first_col[0] == strip_comment_char[0] ) {
            ++first_col;
        }
        if ( col_rename_map.find(first_col) != col_rename_map.end() ) {
            out_line = col_rename_map[first_col].c_str();
        }
        else {
            out_line = first_col;
        }
        if ( strcmp(first_col, colname_feature.c_str()) == 0 ) {
            col_idx_feature = 0;
        }
        if ( strcmp(first_col, colname_x.c_str()) == 0 ) {
            col_idx_x = 0;
        }
        if ( strcmp(first_col, colname_y.c_str()) == 0 ) {
            col_idx_y = 0;
        }
        if ( strcmp(first_col, colname_count.c_str()) == 0 ) {
            col_idx_count = 0;
        }
        for(int32_t i=1; i < tr_mol.nfields; ++i) {
            std::string colname = tr_mol.str_field_at(i);
            // rename the column if specified
            if ( col_rename_map.find(colname) != col_rename_map.end() ) {
                colname = col_rename_map[colname];
            }
            //out_line += (out_mol_tsv_delim + tr_mol.str_field_at(i));
            out_line += (out_mol_tsv_delim + colname);

            if ( strcmp(tr_mol.str_field_at(i), colname_feature.c_str()) == 0 ) {
                if ( col_idx_feature != -1 ) {
                    error("Duplicate column name %s found in the molecule TSV file %s", colname_feature.c_str(), in_mol_tsv.c_str());
                }
                col_idx_feature = i;
            }
            if ( strcmp(tr_mol.str_field_at(i), colname_x.c_str()) == 0 ) {
                if ( col_idx_x != -1 ) {
                    error("Duplicate column name %s found in the molecule TSV file %s", colname_x.c_str(), in_mol_tsv.c_str());
                }
                col_idx_x = i;
            }
            if ( strcmp(tr_mol.str_field_at(i), colname_y.c_str()) == 0 ) {
                if ( col_idx_y != -1 ) {
                    error("Duplicate column name %s found in the molecule TSV file %s", colname_y.c_str(), in_mol_tsv.c_str());
                }
                col_idx_y = i;
            }
            if ( strcmp(tr_mol.str_field_at(i), colname_count.c_str()) == 0 ) {
                if ( col_idx_count != -1 ) {
                    error("Duplicate column name %s found in the molecule TSV file %s", colname_count.c_str(), in_mol_tsv.c_str());
                }
                col_idx_count = i;
            }
        }
        if ( compact_bin ) {
            std::string out_colname_x = col_rename_map.find(colname_x) != col_rename_map.end() ? col_rename_map[colname_x] : colname_x;
            std::string out_colname_y = col_rename_map.find(colname_y) != col_rename_map.end() ? col_rename_map[colname_y] : colname_y;
            std::string out_colname_feature = col_rename_map.find(colname_feature) != col_rename_map.end() ? col_rename_map[colname_feature] : colname_feature;
            std::string out_colname_count = col_rename_map.find(colname_count) != col_rename_map.end() ? col_rename_map[colname_count] : colname_count;
            out_compact_line = out_colname_x + out_mol_tsv_delim + out_colname_y + out_mol_tsv_delim + out_colname_feature + out_mol_tsv_delim + out_colname_count;
        }
    }
    else {
        error("No header line found in the molecule TSV file %s", in_mol_tsv.c_str());
    }

    if ( col_idx_feature == -1 ) {
        error("Column name %s not found in the molecule TSV file %s", colname_feature.c_str(), in_mol_tsv.c_str());
    }
    if ( col_idx_x == -1 ) {
        error("Column name %s not found in the molecule TSV file %s", colname_x.c_str(), in_mol_tsv.c_str());
    }
    if ( col_idx_y == -1 ) {
        error("Column name %s not found in the molecule TSV file %s", colname_y.c_str(), in_mol_tsv.c_str());
    }
    if ( col_idx_count == -1 ) {
        error("Column name %s not found in the molecule TSV file %s", colname_count.c_str(), in_mol_tsv.c_str());
    }

    // print the header line to all output files
    if ( !skip_original ) {
        hprintf(wf_mol_all, "%s\n", out_line.c_str());
    }
    for(int32_t i=0; i < actual_bin_count; ++i) {
        if (compact_bin ) {
            hprintf(wf_mol_bins[i], "%s\n", out_compact_line.c_str());
        }
        else {
            hprintf(wf_mol_bins[i], "%s\n", out_line.c_str());
        }
    }

    uint64_t mol_cnt = 0;
    uint64_t assigned_mol_cnt = 0;
    while ( tr_mol.read_line() > 0 ) {
        std::string ftr_name = tr_mol.str_field_at(col_idx_feature);
        std::map<std::string, int32_t>::iterator it = ftr2bin.find(ftr_name);
        // contruct the output line with renamed columns if specified
        out_line = tr_mol.str_field_at(0);    
        for(int32_t i=1; i < tr_mol.nfields; ++i) {
            out_line += (out_mol_tsv_delim + tr_mol.str_field_at(i));
        }
        if ( !skip_original ) {
            hprintf(wf_mol_all, "%s\n", out_line.c_str());
        }
        if ( it != ftr2bin.end() ) {
            int32_t bin_idx = it->second;
            if ( compact_bin ) {
                hprintf(wf_mol_bins[bin_idx], "%s%s%s%s%s%s%s\n", tr_mol.str_field_at(col_idx_x), out_mol_tsv_delim.c_str(), tr_mol.str_field_at(col_idx_y), out_mol_tsv_delim.c_str(), tr_mol.str_field_at(col_idx_feature), out_mol_tsv_delim.c_str(), tr_mol.str_field_at(col_idx_count));
            }
            else {
                hprintf(wf_mol_bins[bin_idx], "%s\n", out_line.c_str());
            }
            ++assigned_mol_cnt;
        }
        ++mol_cnt;
        if ( mol_cnt % 1000000 == 0 ) {
            notice("%llu molecules processed: (%.2f%%) assigned to bins", mol_cnt, (double)assigned_mol_cnt / mol_cnt * 100);
        }
    }
    if ( assigned_mol_cnt == 0 ) {
        error("No molecules were assigned to any bins. Please check the input files and parameters.");
    }
    notice("%llu molecules processed in total, %llu (%.2f%%) assigned to bins", mol_cnt, assigned_mol_cnt, (double)assigned_mol_cnt / mol_cnt * 100);

    if ( !skip_original ) {
        hts_close(wf_mol_all);
    }
    for(int32_t i=0; i < actual_bin_count; ++i) {
        hts_close(wf_mol_bins[i]);
    }

    notice("Finished writing molecule-level files for each bin");

    // write the feature-level files for each bin
    htsFile* wf_ftr_all = skip_original ? NULL : hts_open((out_prefix + "_all_" + out_ftr_suffix).c_str(), out_ftr_gz ?  "wz" : "w");
    if ( !skip_original ) {
        hprintf(wf_ftr_all, "%s\t%s\n", colname_feature.c_str(), colname_count.c_str());
    }
    std::vector<htsFile*> wf_ftr_bins(actual_bin_count, NULL);
    for(int32_t i=0; i < actual_bin_count; ++i) {
        wf_ftr_bins[i] = hts_open((out_prefix + "_bin" + std::to_string(i+1) + "_" + out_ftr_suffix).c_str(), out_ftr_gz ?  "wz" : "w");
        hprintf(wf_ftr_bins[i], "%s\t%s\n", colname_feature.c_str(), colname_count.c_str());
    }
    for(std::map<std::string, int32_t>::iterator it = ftr2cnts.begin(); it != ftr2cnts.end(); ++it) {
        std::map<std::string, int32_t>::iterator it2 = ftr2bin.find(it->first);
        if ( !skip_original ) {
            hprintf(wf_ftr_all, "%s\t%d\n", it->first.c_str(), it->second);
        }
        if ( it2 != ftr2bin.end() ) {
            int32_t bin_idx = it2->second;
            hprintf(wf_ftr_bins[bin_idx], "%s\t%d\n", it->first.c_str(), it->second);
        }
    }
    hts_close(wf_ftr_all);
    for(int32_t i=0; i < actual_bin_count; ++i) {
        hts_close(wf_ftr_bins[i]);
    }

    notice("Finished writing feature-level files for each bin");

    // write the index file for each bin
    htsFile* wf_index = hts_open((out_prefix + out_index_suffix).c_str(), out_index_gz ? "wz" : "w");
    // take only the basename of the output files for better readability in the index file
    std::string out_prefix_basename = out_prefix;
    size_t last_slash_pos = out_prefix.find_last_of("/\\");
    if ( last_slash_pos != std::string::npos ) {
        out_prefix_basename = out_prefix.substr(last_slash_pos + 1);
    }
    hprintf(wf_index, "bin_id\tmolecule_count\tfeatures_count\tmolecules_path\tfeatures_path\n");
    if ( !skip_original ) {
        hprintf(wf_index, "all\t%llu\t%zu\t%s_all_%s\t%s_all_%s\n", mol_cnt, ftr2cnts.size(), out_prefix_basename.c_str(), out_mol_suffix.c_str(), out_prefix_basename.c_str(), out_ftr_suffix.c_str());
    }
    for(int32_t i=0; i < actual_bin_count; ++i) {
        hprintf(wf_index, "%d\t%llu\t%d\t%s_bin%d_%s\t%s_bin%d_%s\n", i+1, bin_mol_cnts[i], bin_ftr_cnts[i], out_prefix_basename.c_str(), i+1, out_mol_suffix.c_str(), out_prefix_basename.c_str(), i+1, out_ftr_suffix.c_str());
    }
    hts_close(wf_index);
    notice("Finished writing index file");

    // write the JSON file containing bin counts and other metadata
    // Output is a list of "gene", "count", and "bin" for each gene
    std::string out_json_path = out_prefix + out_json_suffix;
    htsFile* wf_json = hts_open(out_json_path.c_str(), "w");
    hprintf(wf_json, "[");
    for(std::map<std::string, int32_t>::iterator it = ftr2cnts.begin();
        it != ftr2cnts.end(); ++it) {
        std::map<std::string, int32_t>::iterator it2 = ftr2bin.find(it->first);
        int32_t bin_idx = it2 != ftr2bin.end() ? it2->second : -1;
        if ( it != ftr2cnts.begin() ) {
            hprintf(wf_json, ",");
        }
        hprintf(wf_json, "{\"gene\":\"%s\",\"count\":%d,\"bin\":%d}", it->first.c_str(), it->second, bin_idx+1);
    }
    hprintf(wf_json, "]\n");
    hts_close(wf_json);
    notice("Finished writing JSON file with bin counts");

    notice("Analysis finished");

    return 0;
}
