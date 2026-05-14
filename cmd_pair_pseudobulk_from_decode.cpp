#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "generic_utils.h"
#include "sge.h"
#include <ctime>
#include <set>

///////////////////////////////////////////////////////////////////////////////////////////////////
// pair-pseudobulk-from-decode : Write a pairwise pseudobulk matrix from pixel-level decode output
///////////////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdPairPseudobulkFromDecode(int32_t argc, char **argv)
{
    std::string tsvf;
    std::string colname_feature("gene");
    std::string colname_count("count");
    std::string colname_factor1;
    std::string colname_factor2;
    std::string missing_value_str = "NA";
    std::string outf;
    std::string suffix_g_f1_f2 = ".g_f1_f2.tsv.gz";
    std::string suffix_f1_f2 = ".f1_f2.tsv.gz";

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("tsv", &tsvf, "tsv file to draw the x-y coordinates. /dev/stdin for stdin")
    LONG_STRING_PARAM("colname-feature", &colname_feature, "Column name for the feature")
    LONG_STRING_PARAM("colname-count", &colname_count, "Column name for the count")
    LONG_STRING_PARAM("colname-factor1", &colname_factor1, "Column name for the first factor")
    LONG_STRING_PARAM("colname-factor2", &colname_factor2, "Column name for the second factor")
    LONG_STRING_PARAM("missing-value-str", &missing_value_str, "Missing value string")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")
    LONG_STRING_PARAM("out-g-f1-f2-suffix", &suffix_g_f1_f2, "Suffix for the output TSV file (default: .g_f1_f2.tsv.gz)")
    LONG_STRING_PARAM("out-f1-f2-suffix", &suffix_f1_f2, "Suffix for the output TSV file (default: .f1_f2.tsv.gz)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( tsvf.empty() || outf.empty() )
        error("--tsv and --out must be specified");

    notice("Analysis started");

    tsv_reader tf(tsvf.c_str());

    // read the header info
    if ( !tf.read_line() )
        error("Cannot read the header from %s", tsvf.c_str());
    std::map<std::string, int32_t> col2idx;
    for(int32_t i=0; i < tf.nfields; ++i) {
        const char* s = tf.str_field_at(i);
        while ( s[0] == '#' ) {
            ++s;
        }
        col2idx[s] = i;
    }

    int32_t icol_feature = find_idx_by_key(col2idx, colname_feature.c_str(), true);
    int32_t icol_count = find_idx_by_key(col2idx, colname_count.c_str(), true);
    int32_t icol_factor1 = find_idx_by_key(col2idx, colname_factor1.c_str(), true);
    int32_t icol_factor2 = find_idx_by_key(col2idx, colname_factor2.c_str(), false);

    if ( icol_feature == -1 ) {
        error("Column name specified by --colname-feature not found in the input TSV file %s", tsvf.c_str());
    }
    if ( icol_count == -1 ) {
        notice("Column name specified by --colname-count not found in the input TSV file %s. Assuming count == 1", tsvf.c_str());
    }
    if ( icol_factor1 == -1 || icol_factor2 == -1 ) {
        error("Both --colname-factor1 and --colname-factor2 must be specified");
    }

    uint64_t nlines = 0;
    std::map<std::string, std::map<std::string, std::map<std::string, double> > > ftr_fac1_fac2_cnt;
    std::map<std::string, std::map<std::string, double> > ftr_fac1_cnt;
    std::map<std::string, std::map<std::string, double> > ftr_fac2_cnt;
    std::map<std::string, std::map<std::string, double> > fac1_fac2_cnt;
    std::map<std::string, double> fac1_cnt;
    std::map<std::string, double> fac2_cnt;
    std::map<std::string, double> ftr_cnt;

    double total = 0;

    while ( tf.read_line() ) {
        const char* ftr = tf.str_field_at(icol_feature);
        int32_t cnt = icol_count >= 0 ? tf.int_field_at(icol_count) : 1;
        const char* fac1 = tf.str_field_at(icol_factor1);
        const char* fac2 = tf.str_field_at(icol_factor2);
        if ( missing_value_str.compare(fac1) == 0 || missing_value_str.compare(fac2) == 0 ) {
            continue;
        }
        ftr_fac1_fac2_cnt[ftr][fac1][fac2] += cnt;
        ftr_cnt[ftr] += cnt;
        fac1_cnt[fac1] += cnt;
        fac2_cnt[fac2] += cnt;
        fac1_fac2_cnt[fac1][fac2] += cnt;
        ftr_fac1_cnt[ftr][fac1] += cnt;
        ftr_fac2_cnt[ftr][fac2] += cnt;
        ftr_fac1_fac2_cnt[ftr][fac1][fac2] += cnt;
        total += cnt;
        ++nlines;
        if ( nlines % 1000000 == 0 ) {
            notice("Processed %d lines", nlines);
        }
    }

    notice("Finished processing %llu lines, identifying %zu factor1, %zu factor2, and across %zu features over %.1f transcripts", nlines, fac1_cnt.size(), fac2_cnt.size(), ftr_cnt.size(), total);

    // write the output count matrix
    std::string out_g_f1_f2 = outf + suffix_g_f1_f2;
    notice("Writing the output count matrix to %s", out_g_f1_f2.c_str());
    htsFile* wf = hts_open(out_g_f1_f2.c_str(), out_g_f1_f2.size() > 3 && out_g_f1_f2.substr(out_g_f1_f2.size() - 3) == ".gz" ? "wz" : "w");
    hprintf(wf, "feature\t%s\t%s\tcount\n", colname_factor1.c_str(), colname_factor2.c_str());
    for( std::map<std::string, std::map<std::string, std::map<std::string, double> > >::iterator it = ftr_fac1_fac2_cnt.begin();
         it != ftr_fac1_fac2_cnt.end(); ++it ) {
        const std::string& ftr = it->first;
        for( std::map<std::string, std::map<std::string, double> >::iterator it2 = it->second.begin();
             it2 != it->second.end(); ++it2 ) {
            const std::string& fac1 = it2->first;
            for( std::map<std::string, double>::iterator it3 = it2->second.begin();
                 it3 != it2->second.end(); ++it3 ) {
                const std::string& fac2 = it3->first;
                double cnt = it3->second;
                if ( cnt == floor(cnt) ) {
                    hprintf(wf, "%s\t%s\t%s\t%.0f\n", ftr.c_str(), fac1.c_str(), fac2.c_str(), cnt);
                }
                else {
                    hprintf(wf, "%s\t%s\t%s\t%.2f\n", ftr.c_str(), fac1.c_str(), fac2.c_str(), cnt);
                }
            }
        }
    }
    hts_close(wf);

    std::string out_f1_f2 = outf + suffix_f1_f2;
    notice("Writing the output factor1-factor2 matrix to %s", out_f1_f2.c_str());
    wf = hts_open(out_f1_f2.c_str(), out_f1_f2.size() > 3 && out_f1_f2.substr(out_f1_f2.size() - 3) == ".gz" ? "wz" : "w");
    hprintf(wf, "%s\t%s\tcount\n", colname_factor1.c_str(), colname_factor2.c_str());
    for( std::map<std::string, std::map<std::string, double> >::iterator it = fac1_fac2_cnt.begin();
         it != fac1_fac2_cnt.end(); ++it ) {
        const std::string& fac1 = it->first;
        for( std::map<std::string, double>::iterator it2 = it->second.begin();
             it2 != it->second.end(); ++it2 ) {
            const std::string& fac2 = it2->first;
            double cnt = it2->second;
            if ( cnt == floor(cnt) ) {
                hprintf(wf, "%s\t%s\t%.0f\n", fac1.c_str(), fac2.c_str(), cnt);
            }
            else {
                hprintf(wf, "%s\t%s\t%.2f\n", fac1.c_str(), fac2.c_str(), cnt);
            }
        }
    }
    hts_close(wf);
    notice("Analysis finished");

    return 0;
}
