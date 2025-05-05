#include "generic_utils.h"
#include "qgenlib/qgen_error.h"
#include <cmath>
#include <cstring>

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
