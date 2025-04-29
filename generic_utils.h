#ifndef __GENERIC_UTILS_H
#define __GENERIC_UTILS_H

#include <string>
#include <map>

// find integer index by key in a dictionary
int32_t find_idx_by_key(std::map<std::string, int32_t>& dict, const char* key, bool required);
double chisq1_log10p(double chi2);

#endif
