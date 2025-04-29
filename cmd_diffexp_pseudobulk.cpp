#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "generic_utils.h"
#include "sge.h"
#include <ctime>
#include <set>

struct deEntry {
    std::string feature;
    std::string factor;
    double chi2;
    double pval;
    double fc;
    double gene_total;
    double log10p;
    deEntry(std::string f, std::string fac, double c, double p, double fch, double g, double l) :
        feature(f), factor(fac), chi2(c), pval(p), fc(fch), gene_total(g), log10p(l) {}
};

bool deEntry_cmp(const deEntry& a, const deEntry& b) {
    return a.log10p > b.log10p;
};

////////////////////////////////////////////////////////////////////////////////
// diffexp-pseudobulk : Pseudobulk test for pixel-level differential expression
////////////////////////////////////////////////////////////////////////////////
int32_t cmdDiffExpPseudobulk(int32_t argc, char **argv)
{
    std::string tsvf;
    std::string colname_feature("feature");
    std::string colname_count("ct");
    std::string colname_K1("K1");
    std::string colname_K2("K2");
    std::string colname_K3("K3");
    std::string colname_P1("P1");
    std::string colname_P2("P2");
    std::string colname_P3("P3");
    double max_pval = 0.001;
    double min_fc = 1.5;
    double pseudocount = 0.5;
    std::string outf;
    std::string suffix_de(".bulk_chisq.tsv");
    std::string suffix_post(".posterior.count.tsv.gz");

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("tsv", &tsvf, "tsv file to draw the x-y coordinates. /dev/stdin for stdin")
    LONG_STRING_PARAM("colname-feature", &colname_feature, "Column name for the feature")
    LONG_STRING_PARAM("colname-count", &colname_count, "Column name for the count")
    LONG_STRING_PARAM("colname-K1", &colname_K1, "Column name for the K1 value")
    LONG_STRING_PARAM("colname-K2", &colname_K2, "Column name for the K2 value")
    LONG_STRING_PARAM("colname-K3", &colname_K3, "Column name for the K3 value")
    LONG_STRING_PARAM("colname-P1", &colname_P1, "Column name for the P1 value")
    LONG_STRING_PARAM("colname-P2", &colname_P2, "Column name for the P2 value")
    LONG_STRING_PARAM("colname-P3", &colname_P3, "Column name for the P3 value")

    LONG_PARAM_GROUP("Settings", NULL)
    LONG_DOUBLE_PARAM("max-pval", &max_pval, "Max p-value for the differential expression test")
    LONG_DOUBLE_PARAM("min-fc", &min_fc, "Min fold change for the differential expression test")
    LONG_DOUBLE_PARAM("pseudocount", &pseudocount, "Pseudocount for the differential expression test")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file name")
    LONG_STRING_PARAM("suffix_de", &suffix_de, "Suffix for differential expression output")
    LONG_STRING_PARAM("suffix_post", &suffix_post, "Suffix for posterior count output")
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
    int32_t icol_K1 = find_idx_by_key(col2idx, colname_K1.c_str(), true);
    int32_t icol_K2 = find_idx_by_key(col2idx, colname_K2.c_str(), false);
    int32_t icol_K3 = find_idx_by_key(col2idx, colname_K3.c_str(), false);
    int32_t icol_P1 = find_idx_by_key(col2idx, colname_P1.c_str(), false);
    int32_t icol_P2 = find_idx_by_key(col2idx, colname_P2.c_str(), false);
    int32_t icol_P3 = find_idx_by_key(col2idx, colname_P3.c_str(), false);

    double min_log10pval = -std::log10(max_pval);
    //double min_log2fc = std::log2(min_fc);

    // possible cases:
    // 1. K1/P1, K2/P2, K3/P3 are present in pairs
    // 2. only K1 present, P1, P2, P3 are empty
    bool best_guess_mode = false;
    int32_t maxK = 0;
    if ( icol_P1 < 0 ) {
        if ( ( icol_P2 >= 0 ) || ( icol_P3 >= 0 ) ) {
            error("When P1 is empty, P2 and P3 must also be empty");
        }
        if ( ( icol_K2 >= 0 ) || ( icol_K3 >= 0 ) ) {
            warning("K2 and K3 is present, but will be ignored");
        }
        best_guess_mode = true;
    }
    else {
        if ( icol_K2 < 0 ) {
            maxK = 1;
        }
        else if ( icol_K3 < 0 ) {
            if ( icol_P2 < 0 ) { error("P2 must present when K2 is present"); }
            maxK = 2;
        }
        else {
            if ( icol_P2 < 0 ) { error("P2 must present when K2 is present"); }
            if ( icol_P3 < 0 ) { error("P3 must present when K3 is present"); }
            maxK = 3;
        }
    }

    uint64_t nlines = 0;
    std::map<std::string, std::map<std::string, double> > fac2ftr2cnt;
    std::map<std::string, double> ftr2cnt;
    std::map<std::string, double> fac2cnt;
    double total = 0;

    std::string out_de = outf + suffix_de;
    std::string out_post = outf + suffix_post;
    htsFile* wf = hts_open(out_de.c_str(), out_de.compare(out_de.size()-3, 3, ".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s", out_de.c_str());
    }
    hprintf(wf, "gene\tfactor\tChi2\tpval\tFoldChange\tgene_total\tlog10p\n");

    while ( tf.read_line() ) {
        const char* ftr = tf.str_field_at(icol_feature);
        int32_t cnt = tf.int_field_at(icol_count);
        const char* K1 = tf.str_field_at(icol_K1);
        if ( best_guess_mode ) {
            fac2ftr2cnt[K1][ftr] += cnt;
            ftr2cnt[ftr] += cnt;
            fac2cnt[K1] += cnt;
            total += cnt;
        }
        else {
            double P1 = tf.double_field_at(icol_P1);
            fac2ftr2cnt[K1][ftr] += (cnt*P1);
            ftr2cnt[ftr] += (cnt*P1);
            fac2cnt[K1] += (cnt*P1);
            total += (cnt*P1);
            if ( maxK >= 2 ) {
                const char* K2 = tf.str_field_at(icol_K2);
                double P2 = tf.double_field_at(icol_P2); 
                fac2ftr2cnt[K2][ftr] += (cnt*P2);
                ftr2cnt[ftr] += (cnt*P2);
                fac2cnt[K2] += (cnt*P2);
                total += (cnt*P2);
                if ( maxK >= 3 ) {
                    const char* K3 = tf.str_field_at(icol_K3);
                    double P3 = tf.double_field_at(icol_P3); 
                    fac2ftr2cnt[K3][ftr] += (cnt*P3);
                    ftr2cnt[ftr] += (cnt*P3);
                    fac2cnt[K3] += (cnt*P3);
                    total += (cnt*P3);
                }
            }
        }
        ++nlines;
        if ( nlines % 1000000 == 0 ) {
            notice("Processed %d lines", nlines);
        }
    }

    notice("Finished processing %llu lines, identifying %zu factors across %zu features over %.1f transcripts", nlines, fac2ftr2cnt.size(), ftr2cnt.size(), total);

    // perform per gene, per factor test
    std::map<std::string, std::map<std::string, double> >::iterator it1;
    std::map<std::string, double>::iterator it2;
    for(it1 = fac2ftr2cnt.begin(); it1 != fac2ftr2cnt.end(); ++it1) {
        std::string fac = it1->first;

        std::vector<deEntry> deEntries;

        for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
            std::string ftr = it2->first;
            double a = it2->second + 0.5;
            double b = fac2cnt[fac] - a + 0.5;
            double c = ftr2cnt[ftr] - a + 0.5;
            double d = total - a - b - c + 0.5;
            double n = a + b + c + d;
            double chisq = (a*d - b*c) * (a*d - b*c) / ((a+b)*(c+d)*(a+c)*(b+d)) * n;
            double log10pval = chisq1_log10p(chisq);
            double fc = (a + pseudocount)*(d + pseudocount) / (b + pseudocount) / (c + pseudocount);

            if ( ( log10pval > min_log10pval ) && ( fc > min_fc ) ) {
                deEntries.push_back(deEntry(ftr, fac, chisq, std::pow(0.1, log10pval), fc, it2->second, log10pval));
            }
         }

        std::sort(deEntries.begin(), deEntries.end(), deEntry_cmp);
        for (size_t i = 0; i < deEntries.size(); ++i) {
            hprintf(wf, "%s\t%s\t%.2f\t%.2e\t%.2f\t%.1f\t%.2f\n",
                    deEntries[i].feature.c_str(), deEntries[i].factor.c_str(),
                    deEntries[i].chi2, deEntries[i].pval, deEntries[i].fc,
                    deEntries[i].gene_total, deEntries[i].log10p);
        }
    }
    hts_close(wf);

    // reverse sort the genes by the count
    std::vector<std::pair<std::string, double> > ftr2cnt_vec;
    for(it2 = ftr2cnt.begin(); it2 != ftr2cnt.end(); ++it2) {
        ftr2cnt_vec.push_back(std::make_pair(it2->first, it2->second));
    }
    std::sort(ftr2cnt_vec.begin(), ftr2cnt_vec.end(), [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
        return a.second > b.second;
    });

    // sort the factors as numbers, and alphanumerically if not numbers
    std::vector<std::string> facs;
    for(it1 = fac2ftr2cnt.begin(); it1 != fac2ftr2cnt.end(); ++it1) {
        facs.push_back(it1->first);
    }
    std::sort(facs.begin(), facs.end(), [](const std::string& a, const std::string& b) {
        if ( a.find_first_not_of("0123456789") == std::string::npos && b.find_first_not_of("0123456789") == std::string::npos ) {
            return std::stoi(a) < std::stoi(b);
        }
        else {
            return a < b;
        }
    });

    // write the posterior count
    notice("Writing posterior count to %s", out_post.c_str());
    htsFile* wf_post = hts_open(out_post.c_str(), out_post.compare(out_post.size()-3, 3, ".gz") == 0 ? "wz" : "w");
    if ( wf_post == NULL ) {
        error("Cannot open output file %s", out_post.c_str());
    }
    hprintf(wf_post, "gene");
    for (size_t i = 0; i < facs.size(); ++i) {
        hprintf(wf_post, "\t%s", facs[i].c_str());
    }
    hprintf(wf_post, "\n");
    for (size_t i = 0; i < ftr2cnt_vec.size(); ++i) {
        std::string ftr = ftr2cnt_vec[i].first;
        hprintf(wf_post, "%s", ftr.c_str());
        for (size_t j = 0; j < facs.size(); ++j) {
            std::string fac = facs[j];
            double cnt = fac2ftr2cnt[fac][ftr];
            hprintf(wf_post, "\t%.2f", cnt);
        }
        hprintf(wf_post, "\n");
    }
    hts_close(wf_post);

    notice("Analysis finished");

    return 0;
}
