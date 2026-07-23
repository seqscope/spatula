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
#include "generic_utils.h"

///////////////////////////////////////////////////////////////////////////
// sptsv-null-residuals : Compute null residuals for a sparse TSV
//
// Under the null model of no gene-by-unit interaction, the expected count
// of gene g in unit u is
//
//   expected_{g,u} = marginal_{g} * marginal_{u} / total
//
// where marginal_{g} is the total count of gene g across all units,
// marginal_{u} is the total count of unit u across all genes, and total is
// the grand total of all counts. This tool reports the L1 deviation from
// that expectation,
//
//   residual_{u} = sum_g |x_{g,u} - expected_{g,u}|   (per unit)
//   residual_{g} = sum_u |x_{g,u} - expected_{g,u}|   (per feature)
//
// Note that the sums run over ALL genes/units, including the (vast majority
// of) zero entries that are not stored in the sparse TSV. Because
// sum_g expected_{g,u} = marginal_{u} and sum_u expected_{g,u} = marginal_{g},
// the contribution of the zero entries can be obtained in closed form from
// the observed entries only:
//
//   residual_{u} = sum_{g in obs(u)} ( |x_{g,u} - E_{g,u}| - E_{g,u} ) + marginal_{u}
//   residual_{g} = sum_{u in obs(g)} ( |x_{g,u} - E_{g,u}| - E_{g,u} ) + marginal_{g}
//
// so a single streaming pass over the sparse TSV suffices.
//////////////////////////////////////////////////////////////////////////
int32_t cmdSptsvNullResiduals(int32_t argc, char **argv)
{
    std::string in_sptsv;   // input sparse TSV prefix
    std::string out_prefix; // output prefix for the residual files
    std::string colname_random_key("random_key");
    int32_t min_feature_count = 0;  // minimum feature count to include in the output
    std::string keyname_dictionary("dictionary");
    std::string keyname_header_info("header_info");
    std::string keyname_n_features("n_features");
    std::string keyname_n_units("n_units");
    std::string keyname_offset_data("offset_data");
    bool exclude_random_key = false;
    std::string suffix_feature_residuals(".feature_residuals_null.tsv");
    std::string suffix_unit_stats(".unit_stats_null.tsv");
    std::string suffix_sptsv_json(".json");
    std::string suffix_feature_stats(".feature.stats.tsv");
    std::string suffix_sptsv_tsv(".txt");

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("sptsv", &in_sptsv, "Input sparse TSV prefix")
    LONG_STRING_PARAM("out", &out_prefix, "Output prefix")
    LONG_INT_PARAM("min-feature-count", &min_feature_count, "Minimum total count for a feature to be included in the null model")

    LONG_PARAM_GROUP("Auxiliary Input/Output Options", NULL)
    LONG_STRING_PARAM("colname-random-key", &colname_random_key, "Column name for the random key in the output")
    LONG_STRING_PARAM("keyname-dictionary", &keyname_dictionary, "Key name for the dictionary in the metadata file")
    LONG_STRING_PARAM("keyname-header-info", &keyname_header_info, "Key name for the header information in the metadata file")
    LONG_STRING_PARAM("keyname-n-features", &keyname_n_features, "Key name for the number of features in the metadata file")
    LONG_STRING_PARAM("keyname-n-units", &keyname_n_units, "Key name for the number of units in the metadata file")
    LONG_STRING_PARAM("keyname-offset-data", &keyname_offset_data, "Key name for the offset data in the metadata file")
    LONG_PARAM("exclude-random-key", &exclude_random_key, "Exclude the random key column in the per-unit output")
    LONG_STRING_PARAM("suffix-feature-residuals", &suffix_feature_residuals, "Suffix for the output per-feature residual file")
    LONG_STRING_PARAM("suffix-unit-stats", &suffix_unit_stats, "Suffix for the output per-unit residual file")
    LONG_STRING_PARAM("suffix-sptsv-json", &suffix_sptsv_json, "Suffix for the input JSON metadata file")
    LONG_STRING_PARAM("suffix-feature-stats", &suffix_feature_stats, "Suffix for the input feature stats file")
    LONG_STRING_PARAM("suffix-sptsv-tsv", &suffix_sptsv_tsv, "Suffix for the input sparse TSV file")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if ( in_sptsv.empty() || out_prefix.empty() ) {
        error("--sptsv and --out must be specified");
    }

    std::string in_ftr_stats = in_sptsv + suffix_feature_stats;
    std::string in_meta = in_sptsv + suffix_sptsv_json;
    std::string in_tsv = in_sptsv + suffix_sptsv_tsv;

    std::map<std::string, int64_t> feature_counts;

    // read the feature stats mapping gene name to feature count
    notice("Reading feature stats file %s", in_ftr_stats.c_str());
    tsv_reader tr_ftr_stats(in_ftr_stats.c_str());
    while ( tr_ftr_stats.read_line() > 0 ) {
        if ( tr_ftr_stats.nfields < 2 ) {
            error("Feature stats file %s has less than 2 columns", in_ftr_stats.c_str());
        }
        const char* feature = tr_ftr_stats.str_field_at(0);
        if ( feature[0] == '#' ) {
            continue; // skip comment lines
        }
        const char* cnt_str = tr_ftr_stats.str_field_at(1);
        if ( ( tr_ftr_stats.nlines == 1 ) && ( strspn(cnt_str, "0123456789") != strlen(cnt_str) ) ) {
            continue; // skip the header line, if any
        }
        feature_counts[feature] = tr_ftr_stats.int64_field_at(1);
    }
    tr_ftr_stats.close();
    notice("Loaded counts for %zu features from %s", feature_counts.size(), in_ftr_stats.c_str());

    // read the metadata
    notice("Reading metadata file %s", in_meta.c_str());
    nlohmann::json json_data;
    int32_t icol_random_key = -1; // column index for the random key, -1 if not present
    int32_t offset_data = -1;     // default offset for the data in the input TSV file
    int32_t n_features = -1;      // number of features in the input TSV file
    int32_t n_units = -1;         // number of units in the input TSV file
    std::vector<std::string> header_info;
    std::vector<std::string> feature_names;
    std::map<std::string, int32_t> ftr2idx;
    {
        std::ifstream meta_file(in_meta);
        if (!meta_file.is_open()) {
            error("Cannot open metadata file %s", in_meta.c_str());
        }
        meta_file >> json_data;

        try {
            header_info = json_data[keyname_header_info].get<std::vector<std::string> >();
            for(int32_t i=0; i < (int32_t)header_info.size(); ++i) {
                if ( header_info[i].compare(colname_random_key) == 0 ) {
                    icol_random_key = i;
                }
            }
        }
        catch (nlohmann::json::exception& e) {
            error("Cannot read header information from the metadata file %s: %s", in_meta.c_str(), e.what());
        }
        offset_data = json_data[keyname_offset_data].get<int32_t>();
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

        if ( icol_random_key < 0 && exclude_random_key ) {
            warning("Random key column %s not found in the metadata file %s. --exclude-random-key has no effect", colname_random_key.c_str(), in_meta.c_str());
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
        if ( (int32_t)header_info.size() != offset_data ) {
            error("The number of header fields (%zu) does not match offset_data (%d) in the metadata file %s", header_info.size(), offset_data, in_meta.c_str());
        }
    }

    notice("Loaded metadata for %d features and %d units from %s", n_features, n_units, in_meta.c_str());

    // set up the marginal counts per feature, applying the feature count threshold
    std::vector<double> ftr_marginals(n_features, 0.0); // marginal_{gene}
    std::vector<bool>   ftr_included(n_features, false);
    double total = 0.0;                                 // total counts across all features/units
    int32_t n_ftr_included = 0;
    int32_t n_ftr_missing = 0;
    for(int32_t i = 0; i < n_features; ++i) {
        std::map<std::string, int64_t>::const_iterator it = feature_counts.find(feature_names[i]);
        if ( it == feature_counts.end() ) {
            ++n_ftr_missing;
            continue; // feature has no count in the feature stats file - always excluded
        }
        ftr_marginals[i] = (double)it->second;
        if ( ( it->second > 0 ) && ( it->second >= (int64_t)min_feature_count ) ) {
            ftr_included[i] = true;
            total += ftr_marginals[i];
            ++n_ftr_included;
        }
    }
    if ( n_ftr_missing > 0 ) {
        warning("%d features in %s are missing from %s, and they will be ignored", n_ftr_missing, in_meta.c_str(), in_ftr_stats.c_str());
    }
    if ( total <= 0 ) {
        error("Total count across the %d included features is zero. Check %s or lower --min-feature-count", n_ftr_included, in_ftr_stats.c_str());
    }
    notice("Using %d of %d features (total counts = %.0lf) to build the null model", n_ftr_included, n_features, total);

    // open the output files
    std::string out_unit_file = out_prefix + suffix_unit_stats;
    std::string out_ftr_file = out_prefix + suffix_feature_residuals;
    htsFile* wf_unit = hts_open(out_unit_file.c_str(), out_unit_file.compare(out_unit_file.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wf_unit == NULL ) {
        error("Cannot open output per-unit residual file %s for writing", out_unit_file.c_str());
    }
    // write the header line, carrying over the header columns of the input sparse TSV
    for(int32_t i = 0; i < offset_data; ++i) {
        if ( ( i == icol_random_key ) && exclude_random_key ) {
            continue;
        }
        hprintf(wf_unit, "%s\t", header_info[i].c_str());
    }
    hprintf(wf_unit, "n_features\tn_umis\tl1_residual\tfrac_residual\n");

    // per-feature accumulators over the observed (non-zero) entries
    std::vector<double>  ftr_sum_abs(n_features, 0.0); // sum_u |x_{g,u} - E_{g,u}|
    std::vector<double>  ftr_sum_exp(n_features, 0.0); // sum_u E_{g,u}
    std::vector<double>  ftr_sum_obs(n_features, 0.0); // sum_u x_{g,u} (to verify the marginals)
    std::vector<int64_t> ftr_n_units(n_features, 0);   // number of units with a non-zero count

    std::vector<int32_t> row_iftrs;
    std::vector<double>  row_counts;
    double sum_unit_marginals = 0.0;
    int64_t n_lines = 0;
    int32_t n_umi_mismatch = 0;

    notice("Reading the sparse TSV file %s", in_tsv.c_str());
    tsv_reader tsv_tr;
    if ( !tsv_tr.open(in_tsv.c_str()) ) {
        error("Cannot open input TSV file %s for reading", in_tsv.c_str());
    }
    while( tsv_tr.read_line() ) {
        int32_t n_genes = tsv_tr.int_field_at(offset_data);
        int64_t n_umis = tsv_tr.int64_field_at(offset_data + 1);
        if ( tsv_tr.nfields != 2 * n_genes + offset_data + 2 ) {
            error("The number of fields in the input TSV file (%d) does not match the expected number (%d) at line %llu",
                  tsv_tr.nfields, 2 * n_genes + offset_data + 2, (unsigned long long)tsv_tr.nlines);
        }

        // first pass over the row : collect the included entries and the unit marginal
        row_iftrs.clear();
        row_counts.clear();
        double unit_marginal = 0.0; // marginal_{unit}
        for(int32_t i = 0; i < n_genes; ++i) {
            int32_t iftr = tsv_tr.int_field_at(2*i + offset_data + 2);
            int32_t cnt  = tsv_tr.int_field_at(2*i + offset_data + 3);
            if ( ( iftr < 0 ) || ( iftr >= n_features ) ) {
                error("Invalid feature index %d in the input TSV file %s at line %llu", iftr, in_tsv.c_str(), (unsigned long long)tsv_tr.nlines);
            }
            if ( !ftr_included[iftr] ) {
                continue;
            }
            row_iftrs.push_back(iftr);
            row_counts.push_back((double)cnt);
            unit_marginal += (double)cnt;
        }

        if ( ( min_feature_count == 0 ) && ( n_ftr_missing == 0 ) && ( (int64_t)unit_marginal != n_umis ) ) {
            ++n_umi_mismatch;
            if ( n_umi_mismatch <= 10 ) {
                warning("The total count of the unit at line %llu (%.0lf) does not match the %dth column of %s (%lld)",
                        (unsigned long long)tsv_tr.nlines, unit_marginal, offset_data + 2, in_tsv.c_str(), (long long)n_umis);
            }
        }
        sum_unit_marginals += unit_marginal;

        // second pass over the row : accumulate the residuals
        // sum_g |x - E| = sum_{g in obs} ( |x - E| - E ) + sum_{all g} E, and sum_{all g} E = marginal_{unit}
        double unit_l1 = unit_marginal;
        for(int32_t i = 0; i < (int32_t)row_iftrs.size(); ++i) {
            int32_t iftr = row_iftrs[i];
            double cnt = row_counts[i];
            double expected = ftr_marginals[iftr] * unit_marginal / total;
            double abs_dev = fabs(cnt - expected);

            unit_l1 += (abs_dev - expected);

            ftr_sum_abs[iftr] += abs_dev;
            ftr_sum_exp[iftr] += expected;
            ftr_sum_obs[iftr] += cnt;
            ++ftr_n_units[iftr];
        }
        if ( unit_l1 < 0 ) unit_l1 = 0; // guard against numerical underflow

        // write the per-unit residuals, carrying over the header columns of the input
        for(int32_t i = 0; i < offset_data; ++i) {
            if ( ( i == icol_random_key ) && exclude_random_key ) {
                continue;
            }
            hprintf(wf_unit, "%s\t", tsv_tr.str_field_at(i));
        }
        // the L1 residual is bounded by 2 * marginal_{unit}, so frac_residual is in [0, 1]
        hprintf(wf_unit, "%zu\t%.0lf\t%.3lf\t%.5lf\n", row_iftrs.size(), unit_marginal, unit_l1,
                unit_marginal > 0 ? unit_l1 / (2.0 * unit_marginal) : 0.0);

        if ( ++n_lines % 1000000 == 0 ) {
            notice("Processed %lld units from %s", (long long)n_lines, in_tsv.c_str());
        }
    }
    tsv_tr.close();
    hts_close(wf_unit);

    notice("Finished processing %lld units, writing %s", (long long)n_lines, out_unit_file.c_str());

    if ( n_umi_mismatch > 0 ) {
        warning("The total count did not match the %dth column of %s for %d units", offset_data + 2, in_tsv.c_str(), n_umi_mismatch);
    }
    if ( fabs(sum_unit_marginals - total) > 1e-6 * total ) {
        warning("The sum of unit marginals (%.0lf) does not match the sum of feature marginals (%.0lf). The null expectation may be inconsistent with %s",
                sum_unit_marginals, total, in_ftr_stats.c_str());
    }

    // write the per-feature residuals
    notice("Writing the per-feature residuals to %s", out_ftr_file.c_str());
    htsFile* wf_ftr = hts_open(out_ftr_file.c_str(), out_ftr_file.compare(out_ftr_file.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wf_ftr == NULL ) {
        error("Cannot open output per-feature residual file %s for writing", out_ftr_file.c_str());
    }
    hprintf(wf_ftr, "feature\tfeature_idx\tn_units\tn_umis\tl1_residual\tfrac_residual\n");
    int32_t n_ftr_mismatch = 0;
    for(int32_t i = 0; i < n_features; ++i) {
        if ( !ftr_included[i] ) {
            continue;
        }
        if ( fabs(ftr_sum_obs[i] - ftr_marginals[i]) > 1e-6 * ftr_marginals[i] ) {
            ++n_ftr_mismatch;
            if ( n_ftr_mismatch <= 10 ) {
                warning("The observed count of feature %s (%.0lf) does not match its count in %s (%.0lf)",
                        feature_names[i].c_str(), ftr_sum_obs[i], in_ftr_stats.c_str(), ftr_marginals[i]);
            }
        }
        // sum_u |x - E| = sum_{u in obs} ( |x - E| - E ) + sum_{all u} E, and sum_{all u} E = marginal_{gene}
        double ftr_l1 = ftr_sum_abs[i] - ftr_sum_exp[i] + ftr_marginals[i];
        if ( ftr_l1 < 0 ) ftr_l1 = 0; // guard against numerical underflow
        hprintf(wf_ftr, "%s\t%d\t%lld\t%.0lf\t%.3lf\t%.5lf\n", feature_names[i].c_str(), i,
                (long long)ftr_n_units[i], ftr_marginals[i], ftr_l1, ftr_l1 / (2.0 * ftr_marginals[i]));
    }
    hts_close(wf_ftr);
    if ( n_ftr_mismatch > 0 ) {
        warning("The observed counts did not match %s for %d features", in_ftr_stats.c_str(), n_ftr_mismatch);
    }

    notice("Analysis finished");

    return 0;
}
