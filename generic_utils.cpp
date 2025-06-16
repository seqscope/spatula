#include "generic_utils.h"
#include <cmath>
#include <cstring>
#include "qgenlib/tsv_reader.h"
#include "qgenlib/hts_utils.h"

int32_t find_idx_by_key(std::map<std::string, int32_t>& dict, const char* key, bool required) {
    auto feature_idx = dict.find(key);
    if ( strlen(key) == 0 ) {
        if ( required ) {
            error("Feature column '%s' cannot be empty", key);
        }
        else {
            return -1;
        }
    }
    if (feature_idx == dict.end()) {
        error("Feature column '%s' not found in the input file", key);
    }
    return feature_idx->second;
}

// compute exp(loga-logb) in a numerically stable way, assuming loga >= logb
double log_sub(double loga, double logb) {
    if ( loga == -std::numeric_limits<double>::infinity() ) return -std::numeric_limits<double>::infinity();
    if ( logb == -std::numeric_limits<double>::infinity() ) return loga;
    if ( loga < logb ) {
        error("loga must be greater than or equal to logb");
    }
    return loga + std::log1p(-std::exp(logb - loga));
}

double log_add(double loga, double logb) {
    if ( loga == -std::numeric_limits<double>::infinity() ) return logb;
    if ( logb == -std::numeric_limits<double>::infinity() ) return loga;
    if ( loga > logb ) {
        return loga + std::log1p(std::exp(logb - loga));
    } else {
        return logb + std::log1p(std::exp(loga - logb));
    }
}

double stable_cauchy_combination_test(const std::vector<double>& log10pvals) {
    std::vector<double> weights(log10pvals.size(), 1.0); // default weights are all 1.0
    return stable_cauchy_combination_test(log10pvals, weights);
}


// returns the CCT log10pval
double stable_cauchy_combination_test(const std::vector<double>& log10pvals, const std::vector<double>& weights) {
    if (log10pvals.size() != weights.size()) {
        error("log10pvals and weights must have the same size");
    }    

    double pi = std::atan(1) * 4;
    double ln_pi = std::log(pi);
    double log10_pi = std::log10(pi);
    double ln_10 = std::log(10.0);
    double log10p_approx_thres = 20.0;

    // normalize the weights
    double sum_weights = 0.0;
    for(int32_t i = 0; i < (int32_t)weights.size(); ++i) {
        sum_weights += weights[i];
    }

    // compute the T statistic, use logT above the threshold, T below the threshold
    double logT = -std::numeric_limits<double>::infinity();
    double T = 0.0; 

    for(int32_t i = 0; i < (int32_t)log10pvals.size(); ++i) {
        if ( weights[i] > 0 ) {
            double w = weights[i] / sum_weights; // normalize the weight
            notice("i = %d, log10pval = %.2f, weight = %.3g", i, log10pvals[i], w);
            if (log10pvals[i] < log10p_approx_thres) {                
                double p = std::pow(10.0, -log10pvals[i]);
                if ( p == 1.0 ) p = 0.999999; // avoid exact 1.0
                T += (w / std::tan(p*pi));
            } 
            else {
                // use approximation : log10T = log10p - log10(pi), 
                // logT = log10p * ln10 - log(pi)
                logT = log_add(logT, log(w/pi) + log10pvals[i] * ln_10);
            }
        }
    }

    double cct_log10pval = 0.0;
    if ( logT == -std::numeric_limits<double>::infinity() ) { // logT was not used, just use T
        // use the standard arctan formula
        if ( T > 0 ) {
            cct_log10pval = -std::log10(std::atan(1/T) / pi);
        }
        else {
            cct_log10pval = -std::log10(0.5 - std::atan(T) / pi);
        }
    }
    else {
        logT += log(1.0 + T/exp(logT)); // add the T term to logT
        cct_log10pval = logT/ln_10 + log10_pi;
    }
    notice("Using logT = %.3g and T = %.3g for CCT, returning %.3f", logT, T, cct_log10pval);
    return cct_log10pval;
}

// p-value from chi-squared value
double chisq1_log10p(double chi2) {
    double z = std::sqrt(chi2 / 2.0);
    if (z < 10.0) {
        double p = std::erfc(z);
        return -std::log10(p);
    } else { // Switch to Taylor expansion
        double logp = -z * z - std::log(z * std::sqrt(M_PI));
        return -logp / std::log(10.0);
    }
}

bool unquote_string(std::string& s) {
    if ( s.size() < 2 ) return false;
    if ( s[0] == '"' && s[s.size()-1] == '"' ) {
        s = s.substr(1, s.size()-2);
        return true;
    }
    else if ( s[0] == '\'' && s[s.size()-1] == '\'' ) {
        s = s.substr(1, s.size()-2);
        return true;
    }
    return false;
}

bool pseudobulk_matrix::load(const char* tsvf) {
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

bool pseudobulk_matrix::add(const pseudobulk_matrix& other) {
    // check if the features and factors exactly match
    bool match_factors = true;
    bool match_features = true;

    if ( factors.size() != other.factors.size() ) {
        match_factors = false;
    } else {
        for (size_t i = 0; i < factors.size(); ++i) {
            if ( factors[i] != other.factors[i] ) {
                match_factors = false;
                break;
            }
        }
    }

    if ( features.size() != other.features.size() ) {
        match_features = false;
    } else {
        for (size_t i = 0; i < features.size(); ++i) {
            if ( features[i] != other.features[i] ) {
                match_features = false;
                break;
            }
        }
    }

    if ( match_factors ) {
        if ( match_features ) {  // simple addition
            for(size_t i = 0; i < features.size(); ++i) {
                for (size_t j = 0; j < factors.size(); ++j) {
                    counts[i][j] += other.counts[i][j];
                }
                rowsums[i] += other.rowsums[i];
            }
            for (size_t j = 0; j < factors.size(); ++j) {
                colsums[j] += other.colsums[j];
            }
            total += other.total;
            return true;
        }
        else { // add additional genes as needed
            int32_t n_new_features = 0;
            for (size_t i = 0; i < other.features.size(); ++i) {
                const std::string& feature = other.features[i];
                auto it = feature2idx.find(feature);
                if ( it == feature2idx.end() ) { // new feature
                    features.push_back(feature);
                    feature2idx[feature] = features.size() - 1;
                    counts.push_back(std::vector<double>(factors.size(), 0.0));
                    for(size_t j = 0; j < factors.size(); ++j) {
                        counts.back()[j] = other.counts[i][j];
                    }
                    rowsums.push_back(other.rowsums[i]);
                }
                else { // existing feature
                    size_t idx = it->second;
                    for (size_t j = 0; j < factors.size(); ++j) {
                        counts[idx][j] += other.counts[i][j];
                    }
                    rowsums[idx] += other.rowsums[i];
                }
            }
            for (size_t j = 0; j < factors.size(); ++j) {
                colsums[j] += other.colsums[j];
            }
            total += other.total;
            if ( n_new_features > 0 ) {
                notice("Features do not match exactly. Added %d new features from the other pseudobulk matrix", n_new_features);
            }
            return true;
        }
    }
    else { // factors do not match
        int32_t n_new_factors = 0;
        int32_t n_new_features = 0;
        std::vector<int32_t> idx_factors;
        for(size_t i=0; i < other.factors.size(); ++i) {
            const std::string& factor = other.factors[i];
            auto it = factor2idx.find(factor);
            if ( it == factor2idx.end() ) { // new factor
                factors.push_back(factor);
                factor2idx[factor] = factors.size() - 1;
                colsums.push_back(0.0);
                n_new_factors++;
                idx_factors.push_back((int32_t)factors.size() - 1);
            }
            else {
                idx_factors.push_back(it->second);
            }
        }
        if ( n_new_factors > 0 ) {
            notice("Factors do not match exactly. Added %d new factors from the other pseudobulk matrix", n_new_factors);
            // resize counts to accommodate new factors
            for (size_t i = 0; i < features.size(); ++i) {
                counts[i].resize(factors.size(), 0.0);
            }
            colsums.resize(factors.size(), 0.0);
        }
        // now add the features
        for (size_t i = 0; i < other.features.size(); ++i) {
            const std::string& feature = other.features[i];
            auto it = feature2idx.find(feature);
            if ( it == feature2idx.end() ) { // new feature
                features.push_back(feature);
                feature2idx[feature] = features.size() - 1;
                counts.push_back(std::vector<double>(factors.size(), 0.0));
                n_new_features++;
                for (size_t j = 0; j < other.factors.size(); ++j) {
                    counts.back()[idx_factors[j]] = other.counts[i][j];
                }
                rowsums.push_back(other.rowsums[i]);
            }
            else { // existing feature
                size_t idx = it->second;
                for (size_t j = 0; j < other.factors.size(); ++j) {
                    counts[idx][idx_factors[j]] += other.counts[i][j];
                }
                rowsums[idx] += other.rowsums[i];
            }
        }
        for (size_t j = 0; j < other.factors.size(); ++j) {
            colsums[idx_factors[j]] += other.colsums[j];
        }
        total += other.total;
        if ( n_new_features > 0 ) {
            notice("Features do not match exactly. Added %d new features from the other pseudobulk matrix", n_new_features);
        }
        return true;
    }
}

bool pseudobulk_matrix::write(const char* filename) {
    htsFile* wf = hts_open(filename, "w");
    if ( wf == NULL) {
        error("Cannot open file %s for writing", filename);
        return false;
    }
    hprintf(wf, "%s\t", feature_name.c_str());
    for (const auto& factor : factors) {
        hprintf(wf, "%s\t", factor.c_str());
    }
    hprintf(wf, "\n");

    for (size_t i = 0; i < features.size(); ++i) {
        hprintf(wf, "%s\t", features[i].c_str());
        for (size_t j = 0; j < factors.size(); ++j) {
            hprintf(wf, "%.4f\t", counts[i][j]);
        }
        hprintf(wf, "\n");
    }
    hts_close(wf);
    notice("Wrote pseudobulk matrix to %s with %zu features and %zu factors, total counts: %.1f", filename, features.size(), factors.size(), total);
    return true;
}
