#ifndef __GENERIC_UTILS_H
#define __GENERIC_UTILS_H

#include <string>
#include <map>
#include "qgenlib/qgen_error.h"

// find integer index by key in a dictionary
int32_t find_idx_by_key(std::map<std::string, int32_t>& dict, const char* key, bool required);
double chisq1_log10p(double chi2);
bool unquote_string(std::string& s);
double stable_cauchy_combination_test(const std::vector<double>& log10pvals, const std::vector<double>& weights);
double stable_cauchy_combination_test(const std::vector<double>& log10pvals);
double log_add(double loga, double logb);
double log_sub(double loga, double logb);

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

    bool load(const char* tsvf); 
    bool add(const pseudobulk_matrix& other);
    bool write(const char* filename);

    pseudobulk_matrix() : total(0.0) {}
    pseudobulk_matrix(const char* tsvf) : total(0.0) {
        if ( !load(tsvf) ) {
            error("Cannot load pseudobulk matrix from %s", tsvf);
        }
    }
};

#endif
