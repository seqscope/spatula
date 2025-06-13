#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "generic_utils.h"
#include "sge.h"
#include <ctime>
#include <map>
#include <vector>

class pseudobulk_matrix {
public:
    std::vector<std::string> features; // features (genes)
    std::vector<std::string> factors;  // factors (K1, K2, K3)
    std::vector< std::vector<double> > counts; // counts for each feature x factor pair
    std::vector<double> rowsums; // row sums (total counts for each feature)
    std::vector<double> colsums; // column sums (total counts for each factor)
    std::map<std::string, int32_t> feature2idx;
    std::map<std::string, int32_t> factor2idx;
    std::string feature_name;
    double total; // total counts across all features and factors

    bool load(const char* tsvf) {
        tsv_reader tf(tsvf);
        if ( !tf.read_line() ) {
            error("Cannot read the header from %s", tsvf);
            return false;
        }

        feature_name.assign(tf.str_field_at(0));
        for(int32_t i=1; i < tf.nfields; ++i) {
            std::string s(tf.str_field_at(i));
            factors.push_back(s);
            factor2idx[s] = i - 1; // store the index of the factor
        }
        colsums.resize(factors.size(), 0.0); // initialize column sums to zero

        total = 0.0;
        while ( tf.read_line() ) {
            double rowsum = 0.0;
            std::string rowname(tf.str_field_at(0));
            features.push_back(rowname);
            feature2idx[rowname] = features.size() - 1; // store the index of the feature
            counts.resize(features.size());
            counts.back().resize(factors.size(), 0.0); // +1 for the feature column
            for(size_t i = 1; i < tf.nfields; ++i) {
                if ( tf.nfields != (int32_t)factors.size() + 1 ) {
                    error("Number of fields in the line (%d) does not match the header (%zu)", tf.nfields, factors.size() + 1);
                    return false;
                }
                double cnt = tf.double_field_at(i);
                counts.back()[i-1] = cnt; // store the count for the feature x factor pair
                rowsum += cnt;
                colsums[i-1] += cnt;
                total += cnt;
            }
            rowsums.push_back(rowsum);
        }
        notice("Loaded pseudobulk matrix from %s with %zu features and %zu factors, total counts: %.1f", tsvf, features.size(), factors.size(), total);
        return true;
    }

    pseudobulk_matrix() : total(0.0) {}
    pseudobulk_matrix(const char* tsvf) : total(0.0) {
        if ( !load(tsvf) ) {
            error("Cannot load pseudobulk matrix from %s", tsvf);
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// diffexp-mode-matrix : Perform differential expression test on the model matrix
////////////////////////////////////////////////////////////////////////////////
int32_t cmdDiffExpModelMatrix(int32_t argc, char **argv)
{
    std::string tsv1f;
    std::string tsv2f;
    std::string outprefix;
    double max_pval = 0.001;
    double min_fc = 1.5;
    double min_count = 10.0;
    double pseudocount = 0.5;
    bool test_pairwise = false;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input/Output options", NULL)
    LONG_STRING_PARAM("tsv1", &tsv1f, "First tsv file containing the model or pseudobulk nmatrix")
    LONG_STRING_PARAM("tsv2", &tsv2f, "(Optional) Second tsv file containing the model or pseudobulk nmatrix")
    LONG_STRING_PARAM("out", &outprefix, "Output file prefix")

    LONG_PARAM_GROUP("Settings", NULL)
    LONG_DOUBLE_PARAM("max-pval", &max_pval, "Max p-value for the differential expression test")
    LONG_DOUBLE_PARAM("min-fc", &min_fc, "Min fold change for the differential expression test")
    LONG_DOUBLE_PARAM("min-count", &min_count, "Minimum observed count for the feature to be considered")

    LONG_DOUBLE_PARAM("pseudocount", &pseudocount, "Pseudocount to add to the counts")
    LONG_PARAM("test-pairwise", &test_pairwise, "Perform pairwise test (1 sample test only)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( tsv1f.empty() || outprefix.empty() )
        error("--tsv1 and --out must be specified");

    notice("Analysis started");

    pseudobulk_matrix pbm1(tsv1f.c_str());
    double log10_max_pval = -std::log10(max_pval);

    if ( tsv2f.empty() ) {
        // perform DE test for each gene x factor pair
        notice("Performing DE test for each gene x factor pair in %s", tsv1f.c_str());
        std::string suffix_marginal(".de.marginal.tsv.gz");
        std::string suffix_pairwise(".de.pairwise.tsv.gz");

        // perform marginal test
        htsFile* wf = hts_open((outprefix + suffix_marginal).c_str(), "wz");
        if ( wf == NULL ) {
            error("Cannot open output file %s", (outprefix + suffix_marginal).c_str());
        }
        hprintf(wf, "%s\tfactor\tChi2\tpval\tFoldChange\tgene_total\tlog10p\n", pbm1.feature_name.c_str());
        for(int32_t i=0; i < (int32_t)pbm1.features.size(); ++i) {
            const std::string& feature = pbm1.features[i];
            for(int32_t j=0; j < (int32_t)pbm1.factors.size(); ++j) {
                double a = pbm1.counts[i][j]; // add pseudocount to avoid zero counts
                double b = pbm1.colsums[j] - a; // total counts for the factor minus the current feature
                double c = pbm1.rowsums[i] - a; // total counts for the feature minus the current factor
                double d = pbm1.total - b - c - a; // total counts minus the counts for the feature and factor

                if ( a < min_count ) {
                    continue; // skip features with low counts
                }

                a += pseudocount;
                b += pseudocount;
                c += pseudocount;
                d += pseudocount;
                double fc = (a * d) / ( b * c );
                if ( fc < min_fc ) {
                    continue; // skip features with low fold change
                }
                double chi2 = (a * d - b * c) * (a * d - b * c) / ((a + b) * (c + d) * (a + c) * (b + d)) * (a + b + c + d);
                double log10pval = chisq1_log10p(chi2);
                if ( log10pval < log10_max_pval ) {
                    continue; // skip features with high p-value
                }
                const std::string& factor = pbm1.factors[j];
                // print the results
                hprintf(wf, "%s\t%s\t%.2f\t%.2e\t%.2f\t%.1f\t%.2f\n",
                        feature.c_str(), factor.c_str(),
                        chi2, std::pow(10.0, -log10pval), fc,
                        pbm1.rowsums[i], log10pval);
            }
        }
        hts_close(wf);

        notice("Marginal DE test completed, results written to %s", (outprefix + suffix_marginal).c_str());

        // perform pairwise test
        if ( test_pairwise ) {
            notice("Peforming pairwise DE test for each pairs of factors in %s", tsv1f.c_str());
            htsFile* wf = hts_open((outprefix + suffix_pairwise).c_str(), "wz");
            if ( wf == NULL ) {
                error("Cannot open output file %s", (outprefix + suffix_marginal).c_str());
            }
            hprintf(wf, "%s\tfactor1\tfactor2\tChi2\tpval\tFoldChange\tgene_total\tlog10p\tfrac1\tfrac2\n", pbm1.feature_name.c_str());
            for(int32_t i=0; i < (int32_t)pbm1.features.size(); ++i) {
                const std::string& feature = pbm1.features[i];
                for(int32_t j=1; j < (int32_t)pbm1.factors.size(); ++j) {
                    const std::string& factor1 = pbm1.factors[j];
                    for(int32_t k=0; k < j; ++k) {
                        double a = pbm1.counts[i][j];
                        double b = pbm1.counts[i][k];
                        double c = pbm1.colsums[j] - a;
                        double d = pbm1.colsums[k] - b;

                        if ( a < min_count || b < min_count ) {
                            continue; // skip features with low counts
                        }

                        a += pseudocount; // add pseudocount to avoid zero counts
                        b += pseudocount;
                        c += pseudocount;
                        d += pseudocount;

                        double fc = (a * d) / (b * c);
                        bool swapped = false;
                        if ( fc < 1.0 ) {
                            double tmp = a;
                            a = b;
                            b = tmp;
                            tmp = c;
                            c = d;
                            d = tmp;
                            fc = 1.0 / fc;
                            swapped = true;
                        }
                        if ( fc < min_fc ) {
                            continue; // skip features with low fold change
                        }
                        a += pseudocount;
                        b += pseudocount;
                        c += pseudocount;
                        d += pseudocount;
                        double chi2 = (a * d - b * c) * (a * d - b * c) / ((a + b) * (c + d) * (a + c) * (b + d)) * (a + b + c + d);
                        double log10pval = chisq1_log10p(chi2);
                        if ( log10pval < log10_max_pval ) {
                            continue; // skip features with high p-value
                        }
                        // print the results
                        const std::string& factor2 = pbm1.factors[k];
                        hprintf(wf, "%s\t%s\t%s\t%.2f\t%.2e\t%.2f\t%.1f\t%.2f\t%.5g\t%.5g\n",
                                feature.c_str(), 
                                swapped ? factor2.c_str() : factor1.c_str(),
                                swapped ? factor1.c_str() : factor2.c_str(),
                                chi2, std::pow(10.0, -log10pval), fc,
                                a + b - pseudocount * 2, log10pval, a/(a+c), b/(b+d));
                    }
                }
            }
            hts_close(wf);
        }
    }
    else {
        error("Tests with two pseudobulk matrices are not implemented yet");
        //pseudobulk_matrix pbm2(tsv2f.c_str());
    }

    notice("Analysis finished");

    return 0;
}
