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
    bool ignore_mismatch = false;

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
    LONG_PARAM("ignore-mismatch", &ignore_mismatch, "Ignore mismatching factors between tsv1 and tsv2. If set, only overlapping factors will be used for the test")
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
        //error("Tests with two pseudobulk matrices are not implemented yet");
        pseudobulk_matrix pbm2(tsv2f.c_str());

        // make sure that the two matrices are compatible, i.e. they have the same features and factors
        if ( pbm1.features.size() != pbm2.features.size() || pbm1.factors.size() != pbm2.factors.size() ) {
            if ( ignore_mismatch) {
                notice("WARNING: The two pseudobulk matrices have different number of features or factors, but --ignore-mismatch is set. Only overlapping features and factors will be used for the test.");
            }
            else {
                error("The two pseudobulk matrices must have the same number of features and factors");
            }
        }

        bool match_features = true;
        bool match_factors = true;
        for (size_t i = 0; i < pbm1.features.size(); ++i) {
            if ( pbm1.features[i] != pbm2.features[i] ) {
                match_features = false;
                break;
            }
        }
        for (size_t i = 0; i < pbm1.factors.size(); ++i) {
            if ( pbm1.factors[i] != pbm2.factors[i] ) {
                match_factors = false;
                break;
            }
        }

        // order of index to match tsv2 to tsv1
        std::vector<int32_t> idx1_factors(pbm1.factors.size(), -1);
//        std::vector<int32_t> idx2_factors(pbm2.factors.size(), -1);
        std::vector<int32_t> idx2_factors(pbm1.factors.size(), -1);
        if ( match_factors ) {
            for(size_t i = 0; i < pbm1.factors.size(); ++i) {
                idx1_factors[i] = i;
                idx2_factors[i] = i;
            }
        }
        else {
            for(size_t i = 0; i < pbm1.factors.size(); ++i) {
                idx1_factors[i] = i;
                auto it = pbm2.factor2idx.find(pbm1.factors[i]);
                if ( it == pbm2.factor2idx.end() ) {
                    if ( ignore_mismatch ) {
                        idx2_factors[i] = -1; // skip factors that are not present in the second matrix
                    }
                    else {
                        error("Factor %s from the first pseudobulk matrix not found in the second pseudobulk matrix", pbm1.factors[i].c_str());
                    }
                }
                else {
                    idx2_factors[i] = it->second;
                }
            }
        }
        std::vector<int32_t> idx1_features(pbm1.features.size(), -1);
//        std::vector<int32_t> idx2_features(pbm2.features.size(), -1);
        std::vector<int32_t> idx2_features(pbm1.features.size(), -1);
        if ( match_features ) {
            for(size_t i = 0; i < pbm1.features.size(); ++i) {
                idx1_features[i] = i;
                idx2_features[i] = i;
            }
        }
        else {
            for(size_t i = 0; i < pbm1.features.size(); ++i) {
                idx1_features[i] = i;
                auto it = pbm2.feature2idx.find(pbm1.features[i]);
                if ( it == pbm2.feature2idx.end() ) {
                    if ( ignore_mismatch ) {
                        idx2_features[i] = -1; // skip features that are not present in the second matrix
                    }
                    else {
                        error("Feature %s from the first pseudobulk matrix not found in the second pseudobulk matrix", pbm1.features[i].c_str());
                    } 
                }
                else {
                    idx2_features[i] = it->second;
                }
            }
        }

        // perform the following DE tests:
        // 1. de.bulk.factor marginal test for each factor (colsums)
        // 2. de.bulk.feature marginal test for each feature (rowsums)
        // 3. de.cond.factor conditional test for each factor given a gene
        // 4. de.cond.feature conditional test for each feature given a factor
        // 5. de.cct.factor cauchy-combined test for each factor given a gene
        // 6. de.cct.feature cauchy-combined test for each feature given a factor
        notice("Performing DE tests for each factor and feature");
        
        // de.bulk.factor
        std::string suffix_bulk_factor(".de.bulk.factor.tsv.gz");
        htsFile* wf_bulk_factor = hts_open((outprefix + suffix_bulk_factor).c_str(), "wz");
        if ( wf_bulk_factor == NULL ) {
            error("Cannot open output file %s", (outprefix + suffix_bulk_factor).c_str());
        }
        hprintf(wf_bulk_factor, "Factor\tCount1\tCount2\tFrac1\tFrac2\tlog2FC\tChi2\tpval\tlog10p\n");
        double log2fc_thres = std::log2(min_fc);
        for(int32_t i=0; i < (int32_t)idx1_factors.size(); ++i) {
            int32_t i1 = idx1_factors[i];
            int32_t i2 = idx2_factors[i];
            if ( i1 < 0 || i2 < 0 ) {
                // skip factors that are not present in both matrices
                continue;
            }
            double a = pbm1.colsums[i1];
            double b = pbm2.colsums[i2];
            double c = pbm1.total - a;
            double d = pbm2.total - b;
            if ( a < min_count || b < min_count ) {
                continue; // skip factors with low counts
            }
            a += pseudocount; // add pseudocount to avoid zero counts
            b += pseudocount;
            c += pseudocount;
            d += pseudocount;
            double log2fc = log2((a * d) / (b * c));
            // if ( fabs(log2fc) < log2fc_thres ) {
            //     continue; // skip factors with low fold change
            // }
            double chi2 = (a * d - b * c) * (a * d - b * c) / ((a + b) * (c + d) * (a + c) * (b + d)) * (a + b + c + d);
            double log10pval = chisq1_log10p(chi2);
            // if ( log10pval < log10_max_pval ) {
            //     continue; // skip factors with high p-value
            // }
            // print the results
            hprintf(wf_bulk_factor, "%s\t%.2f\t%.2f\t%.5g\t%.5g\t%.4f\t%.2f\t%.2e\t%.2f\n",
                    pbm1.factors[i].c_str(),
                    a, b,
                    a / pbm1.total, b / pbm2.total,
                    log2fc, chi2, std::pow(10.0, -log10pval), log10pval);
        }
        hts_close(wf_bulk_factor);

        // de.bulk.feature
        std::string suffix_bulk_feature(".de.bulk.feature.tsv.gz");
        htsFile* wf_bulk_feature = hts_open((outprefix + suffix_bulk_feature).c_str(), "wz");
        if ( wf_bulk_feature == NULL ) {
            error("Cannot open output file %s", (outprefix + suffix_bulk_feature).c_str());
        }
        hprintf(wf_bulk_feature, "Feature\tCount1\tCount2\tFrac1\tFrac2\tlog2FC\tChi2\tpval\tlog10p\n");
        for(int32_t i=0; i < (int32_t)idx1_features.size(); ++i) {
            int32_t i1 = idx1_features[i];
            int32_t i2 = idx2_features[i];
            if ( i1 < 0 || i2 < 0 ) {
                // skip features that are not present in both matrices
                continue;
            }
            double a = pbm1.rowsums[i1];
            double b = pbm2.rowsums[i2];
            double c = pbm1.total - a;
            double d = pbm2.total - b;
            if ( a < min_count || b < min_count ) {
                continue; // skip features with low counts
            }
            a += pseudocount; // add pseudocount to avoid zero counts
            b += pseudocount;
            c += pseudocount;
            d += pseudocount;
            double log2fc = log2((a * d) / (b * c));
            if ( fabs(log2fc) < log2fc_thres ) {
                continue; // skip features with low fold change
            }
            double chi2 = (a * d - b * c) * (a * d - b * c) / ((a + b) * (c + d) * (a + c) * (b + d)) * (a + b + c + d);
            double log10pval = chisq1_log10p(chi2);
            if ( log10pval < log10_max_pval ) {
                continue; // skip features with high p-value
            }
            // print the results
            hprintf(wf_bulk_feature, "%s\t%.2f\t%.2f\t%.5g\t%.5g\t%.4f\t%.2f\t%.2e\t%.2f\n",
                    pbm1.features[i].c_str(),
                    a, b,
                    a / pbm1.total, b / pbm2.total,
                    log2fc, chi2, std::pow(10.0, -log10pval), log10pval);
        }
        hts_close(wf_bulk_feature);
        notice("Bulk DE tests completed, results written to %s and %s", 
               (outprefix + suffix_bulk_factor).c_str(), 
               (outprefix + suffix_bulk_feature).c_str());

        // de.cond.feature
        // de.combinedfeature
        std::string suffix_conditional_feature(".de.conditional.feature.tsv.gz");
        std::string suffix_combined_feature(".de.combined.feature.tsv.gz");
        htsFile* wf_conditional_feature = hts_open((outprefix + suffix_conditional_feature).c_str(), "wz");
        if ( wf_conditional_feature == NULL ) {
            error("Cannot open output file %s", (outprefix + suffix_conditional_feature).c_str());
        }
        htsFile* wf_combined_feature = hts_open((outprefix + suffix_combined_feature).c_str(), "wz");
        if ( wf_combined_feature == NULL ) {
            error("Cannot open output file %s", (outprefix + suffix_combined_feature).c_str());
        }
        hprintf(wf_conditional_feature, "Feature\tFactor\tCount1\tCount2\tFrac1\tFrac2\tlog2FC\tChi2\tpval\tlog10p\n");
        hprintf(wf_combined_feature, "Feature\tCount1\tCount2\tFrac1\tFrac2\tTopFactor\tTopLog2FC\tTopLog10p\tCombinedLog10p\n");
        for(int32_t i=0; i < (int32_t)idx1_features.size(); ++i) {
            std::vector<double> log10pvals(idx1_factors.size(), 0.0);
            std::vector<double> weights(idx1_factors.size(), 0.0);
            double top_log10pval = 0;
            double top_log10fc = 0;
            int32_t top_factor_idx = -1;
            for(int32_t j=0; j < (int32_t)idx1_factors.size(); ++j) {
                int32_t i1 = idx1_features[i];
                int32_t j1 = idx1_factors[j];
                int32_t i2 = idx2_features[i];
                int32_t j2 = idx2_factors[j];
                if ( i1 < 0 || j1 < 0 || i2 < 0 || j2 < 0 ) {
                    // skip features or factors that are not present in both matrices
                    continue;
                }
                double a = pbm1.counts[i1][j1];
                double b = pbm2.counts[i2][j2];
                double c = pbm1.colsums[j1] - a;
                double d = pbm2.colsums[j2] - b;

                if ( a < min_count || b < min_count ) {
                    continue; // skip features with low counts
                }

                a += pseudocount; // add pseudocount to avoid zero counts
                b += pseudocount;
                c += pseudocount;
                d += pseudocount;
                double log2fc = log2((a * d) / (b * c));

                double chi2 = (a * d - b * c) * (a * d - b * c) / ((a + b) * (c + d) * (a + c) * (b + d)) * (a + b + c + d);
                double log10pval = chisq1_log10p(chi2);

                log10pvals[j] = log10pval;
                weights[j] = 1.0;

                if ( log10pval > top_log10pval ) {
                    top_log10pval = log10pval;
                    top_log10fc = log2fc;
                    top_factor_idx = j;
                }

                if ( fabs(log2fc) < log2fc_thres ) {
                    continue; // skip features with low fold change
                }
                if ( log10pval < log10_max_pval ) {
                    continue; // skip features with high p-value
                }

                // print the results
                hprintf(wf_conditional_feature, "%s\t%s\t%.2f\t%.2f\t%.5g\t%.5g\t%.4f\t%.2f\t%.2e\t%.2f\n",
                        pbm1.features[i].c_str(), pbm1.factors[j].c_str(),
                        a, b,
                        a / pbm1.total, b / pbm2.total,
                        log2fc, chi2, std::pow(10.0, -log10pval), log10pval);

            }
            // perform cauchy-combined test
            if ( top_factor_idx >= 0 ) {
                notice("Computing CCT for feature %s", pbm1.features[i].c_str());
                double combined_log10pval = stable_cauchy_combination_test(log10pvals, weights);
                // if ( combined_log10pval < log10_max_pval ) {
                //     continue; // skip features with high combined p-value
                // }
                hprintf(wf_combined_feature, "%s\t%.2f\t%.2f\t%.5g\t%.5g\t%s\t%.4f\t%.2f\t%.2f\n",
                        pbm1.features[i].c_str(),
                        pbm1.rowsums[i], pbm2.rowsums[i],
                        pbm1.rowsums[i] / pbm1.total, pbm2.rowsums[i] / pbm2.total,
                        pbm1.factors[top_factor_idx].c_str(),
                        top_log10fc, top_log10pval, combined_log10pval);
            }
        }
        hts_close(wf_conditional_feature);
        hts_close(wf_combined_feature);

        notice("Conditional and Combined DE tests completed, results written to %s and %s", 
               (outprefix + suffix_conditional_feature).c_str(), 
               (outprefix + suffix_combined_feature).c_str());
    }

    notice("Analysis finished");

    return 0;
}
