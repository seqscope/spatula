#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include <cmath>
#include <ctime>
#include <climits>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>

/////////////////////////////////////////////////////////////////////////////////////////
// assign-feature2bin : Assign each feature (gene) into ordered bins based on total counts
//
// This is the first of the two commands split out from the (deprecated)
// 'split-molecule-counts' command. It only performs the gene->bin assignment step,
// writing the resulting assignment as a JSON file (equivalent to the old
// '_bin_counts.json' output). The companion command 'split-mol2bin' consumes this JSON
// to split a molecule-level TSV file into per-bin files.
/////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdAssignFeature2Bin(int32_t argc, char **argv)
{
    std::string in_ftr_tsv;  // TSV file containing gene-level total counts
    std::string out_json;    // Output JSON file containing bin counts and assignment
    std::string in_ftr_tsv_delim = "\t"; // Delimiter for the input feature TSV file
    std::string strip_comment_char = "#"; // Character to strip from the beginning of comment lines
    int32_t bin_count = 50;  // Number of bins to split the data into (roughly equal number of molecules per bin)

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("feature-tsv", &in_ftr_tsv, "TSV file containing gene-level total counts (gene name in column 1, count in column 2)")
    LONG_STRING_PARAM("out-json", &out_json, "Output JSON file containing the gene->bin assignment (equivalent to the old _bin_counts.json)")

    LONG_PARAM_GROUP("Key Parameters", NULL)
    LONG_INT_PARAM("bin-count", &bin_count, "Number of bins to split the features into (roughly equal number of molecules per bin)")

    LONG_PARAM_GROUP("Auxilary Input/Output Parameters", NULL)
    LONG_STRING_PARAM("in-feature-tsv-delim", &in_ftr_tsv_delim, "Delimiter for the input feature TSV file")
    LONG_STRING_PARAM("strip-comment-char", &strip_comment_char, "Character to strip from the beginning of comment lines in the input file (if any)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if ( in_ftr_tsv.empty() ) {
        error("Missing required option --feature-tsv");
    }
    if ( out_json.empty() ) {
        error("Missing required option --out-json");
    }

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

    // write the JSON file containing bin counts and other metadata
    // Output is a list of "gene", "count", and "bin" for each gene
    notice("Writing JSON file with bin counts to %s", out_json.c_str());
    htsFile* wf_json = hts_open(out_json.c_str(), "w");
    if ( wf_json == NULL ) {
        error("Cannot open output JSON file %s for writing", out_json.c_str());
    }
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
