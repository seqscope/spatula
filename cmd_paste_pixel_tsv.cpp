#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge2.h"
#include "file_utils.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

struct _pix_factor_t {
    double x;
    double y;
    std::vector<uint8_t> factors;
    std::vector<double> probs;
    _pix_factor_t(tsv_reader* ptr, int32_t icol_x, int32_t icol_y, int32_t icol_k1, int32_t icol_p1, int32_t topk, double offset_x, double offset_y, double scale) {
        x = ptr->double_field_at(icol_x)/scale + offset_x;
        y = ptr->double_field_at(icol_y)/scale + offset_y;
        for(int32_t i = 0; i < topk; ++i) {
            factors.push_back((uint8_t)ptr->int_field_at(icol_k1 + i));
            probs.push_back(ptr->double_field_at(icol_p1 + i));
        }
    }
};

typedef struct _pix_factor_t pix_factor_t;

/////////////////////////////////////////////////////////////////////////////////////////
// join-pixel-tsv : Join FICTURE's pixel-level output with raw transcript-level TSV files
/////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdPastePixelTSV(int32_t argc, char **argv)
{
    std::string in_mol_tsv;   // TSV file containing individual molecules
    std::vector<std::string> pix_prefix_tsvs; // vector of "[prefix],[tsv-path]" pairs for pixel-level projections
    std::string out_tsv;     // Output TSV file name

    std::string out_suffix_tsv(".tsv.gz");            // suffix for the output TSV file

    std::string colname_X("X");
    std::string colname_Y("Y");
    std::string csv_colnames_include; // comma-separated column names to include in the output TSV file
    std::string csv_colnames_exclude; // comma-separated column names to exclude in the output TSV file

    int32_t out_max_k = 1;            // maximum num of pixel-level factors to include in the joined output 
    int32_t out_max_p = 1;            // maximum num of pixel-level probs to include in the joined output

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_MULTI_STRING_PARAM("pix-prefix-tsv", &pix_prefix_tsvs, "TSV file containing pixel-level factors")
    LONG_STRING_PARAM("out-tsv", &out_tsv, "Output TSV file")

    LONG_PARAM_GROUP("Expected columns in input and output", NULL)
    LONG_STRING_PARAM("colname-x", &colname_X, "Column name for X-axis")
    LONG_STRING_PARAM("colname-y", &colname_Y, "Column name for Y-axis")
    LONG_STRING_PARAM("colnames-include", &csv_colnames_include, "Comma-separated column names to include in the output TSV file")
    LONG_STRING_PARAM("colnames-exclude", &csv_colnames_exclude, "Comma-separated column names to exclude in the output TSV file")
    LONG_INT_PARAM("out-max-k", &out_max_k, "Maximum number of pixel-level factors to include in the joined output. (Default : 1)")
    LONG_INT_PARAM("out-max-p", &out_max_p, "Maximum number of pixel-level posterior probabilities to include in the joined output. (Default : 1)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check the input files
    if ( pix_prefix_tsvs.empty() || out_tsv.empty() )
        error("--pix-prefix-tsv, and --out-prefix must be specified");

    // parse pix_prefix_tsvs
    std::vector<std::string> pix_prefixes;
    std::vector<std::string> pix_tsvs;
    for(int32_t i=0; i < (int32_t)pix_prefix_tsvs.size(); ++i) {
        //notice("Parsing %s..", pix_prefix_tsvs[i].c_str());
        std::vector<std::string> toks;
        split(toks, ",", pix_prefix_tsvs[i].c_str());
        //notice("toks.size = %zu", toks.size());
        if ( toks.size() != 2 ) {
            error("Cannot parse %s", pix_prefix_tsvs[i].c_str());
        }
        pix_prefixes.push_back(toks[0]);
        pix_tsvs.push_back(toks[1]);
    }
    int32_t n_pix = (int32_t)pix_prefixes.size();
    //notice("n_pix = %d", n_pix);

    // read the meta/header lines of the pixel-level factors
    std::vector<tsv_reader*> pix_trs;
    for(int32_t i=0; i < n_pix; ++i) {
        tsv_reader* ptr = new tsv_reader(pix_tsvs[i].c_str());
        pix_trs.push_back(ptr);
    }

    bool mol_colnames_exclude = true;
    std::set<std::string> mol_colnames_set;
    if ( csv_colnames_exclude.empty() ) {
        if ( csv_colnames_include.empty() ) {
            mol_colnames_exclude = true;
        }
        else {
            mol_colnames_exclude = false;
            std::vector<std::string> colnames;
            split(colnames, ",", csv_colnames_include);
            for(int32_t i=0; i < (int32_t)colnames.size(); ++i) {
                mol_colnames_set.insert(colnames[i]);
            }
        }
    }
    else if ( csv_colnames_include.empty() ) {
        mol_colnames_exclude = true;
        std::vector<std::string> colnames;
        split(colnames, ",", csv_colnames_exclude);
        for(int32_t i=0; i < (int32_t)colnames.size(); ++i) {
            mol_colnames_set.insert(colnames[i]);
        }
    }
    else {
        error("Cannot specify both --csv-colnames-exclude and --csv-colnames-include");
    }

    htsFile* wh_tsv = hts_open(out_tsv.c_str(), out_tsv.compare(out_tsv.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");

    bool is_header = true;
    std::vector< int32_t > icols_x;
    std::vector< int32_t > icols_y;
    std::vector< int32_t > icols_include;
    std::vector< std::vector< int32_t > > icols_factor_K_idxs(n_pix);
    std::vector< std::vector< int32_t > > icols_factor_P_idxs(n_pix);
    uint64_t nlines = 0;
    while ( pix_trs[0]->read_line() ) {
        for( int32_t i=1; i < n_pix; ++i ) {
            pix_trs[i]->read_line();
        }
        if ( is_header ) {
            // determine the column indices to include
            for(int32_t i=0; i < n_pix; ++i) {
                std::map<std::string, int32_t> col2idx;
                for(int32_t j=0; j < pix_trs[i]->nfields; ++j) {
                    const char* s = pix_trs[i]->str_field_at(j);
                    if ( s[0] == '#' ) {
                        col2idx[s+1] = j;
                        //notice("col2idx[%s] = %d", s+1, j);
                    }
                    else {
                        col2idx[s] = j;
                        //notice("col2idx[%s] = %d", s, j);
                    }
                }
                if ( col2idx.find(colname_X) == col2idx.end() ) {
                    error("Cannot find %s in the header line of %s", colname_X.c_str(), pix_tsvs[i].c_str());
                }
                icols_x.push_back(col2idx[colname_X]);
                if ( col2idx.find(colname_Y) == col2idx.end() ) {
                    error("Cannot find %s in the header line of %s", colname_Y.c_str(), pix_tsvs[i].c_str());
                }
                icols_y.push_back(col2idx[colname_Y]);
                //notice("i = %d, icols_x = %d, icols_y = %d", i, icols_x.back(), icols_y.back());
                for(int32_t j=0; j < out_max_k; ++j) {
                    std::string colname = "K" + std::to_string(j+1);
                    if ( col2idx.find(colname) == col2idx.end() ) {
                        error("Cannot find %s in the header line of %s", colname.c_str(), pix_tsvs[i].c_str());
                    }
                    icols_factor_K_idxs[i].push_back(col2idx[colname]);
                }
                for(int32_t j=0; j < out_max_p; ++j) {
                    std::string colname = "P" + std::to_string(j+1);
                    if ( col2idx.find(colname) == col2idx.end() ) {
                        error("Cannot find %s in the header line of %s", colname.c_str(), pix_tsvs[i].c_str());
                    }
                    icols_factor_P_idxs[i].push_back(col2idx[colname]);
                }
                if ( i == 0 ) {
                    if ( mol_colnames_exclude ) {
                        mol_colnames_set.insert(colname_X);
                        mol_colnames_set.insert(colname_Y);
                        int32_t icol_K1 = -1;
                        std::string colname("K1");
                        if ( col2idx.find(colname) == col2idx.end() ) {
                            icol_K1 = pix_trs[i]->nfields;
                        }
                        else {
                            icol_K1 = col2idx[colname];
                        }
                        // add the column names to include
                        for(int32_t j=0; j < icol_K1; ++j) {
                            const char* s = pix_trs[i]->str_field_at(j);
                            if ( s[0] == '#' ) ++s;
                            if ( mol_colnames_set.find(s) == mol_colnames_set.end() ) {
                                icols_include.push_back(j);
                            }
                        }
                    }
                    else {
                        for(int32_t j=0; j < pix_trs[i]->nfields; ++j) {
                            const char* s = pix_trs[i]->str_field_at(j);
                            if ( s[0] == '#' ) ++s; 
                            if ( mol_colnames_set.find(s) != mol_colnames_set.end() ) {
                                icols_include.push_back(j);
                            }
                        }
                    }
                    // print out the header lines to include
                    hprintf(wh_tsv, "%s\t%s", colname_X.c_str(), colname_Y.c_str());
                    for(int32_t j=0; j < icols_include.size(); ++j) {
                        hprintf(wh_tsv, "\t%s", pix_trs[i]->str_field_at(icols_include[j]));
                    }
                }
                for(int32_t j=0; j < out_max_k; ++j) {
                    hprintf(wh_tsv, "\t%sK%d", pix_prefixes[i].c_str(), j+1);
                }
                for(int32_t j=0; j < out_max_p; ++j) {
                    hprintf(wh_tsv, "\t%sP%d", pix_prefixes[i].c_str(), j+1);
                }
            }
            hprintf(wh_tsv, "\n");
            is_header = false;
        }
        else {
            // print out the actual data lines
            //notice("icols_x[0] = %d, icols_y[0] = %d", icols_x[0], icols_y[0]);
            double mol_x = 0.0, mol_y = 0.0;
            for(int32_t i=0; i < n_pix; ++i) {
                if ( i == 0 ) {
                    hprintf(wh_tsv, "%s", pix_trs[i]->str_field_at(icols_x[i]));
                    hprintf(wh_tsv, "\t%s", pix_trs[i]->str_field_at(icols_y[i]));
                    for(int32_t j=0; j < icols_include.size(); ++j) {
                        hprintf(wh_tsv, "\t%s", pix_trs[i]->str_field_at(icols_include[j]));
                    }
                    mol_x = pix_trs[i]->double_field_at(icols_x[i]);
                    mol_y = pix_trs[i]->double_field_at(icols_y[i]);
                    if ( icols_x[i] == icols_y[i] ) {
                        error("X and Y coordinates uses the same column %d at i = %d", icols_x[i], i);
                    }
                }
                else {
                    if ( mol_x != pix_trs[i]->double_field_at(icols_x[i]) ) {
                        error("(%f, %f) coordinates are different from (%.f, %f) at line %llu", pix_trs[i]->double_field_at(icols_x[i]), pix_trs[i]->double_field_at(icols_y[i]), mol_x, mol_y, nlines);
                    }
                    if ( mol_y != pix_trs[i]->double_field_at(icols_y[i]) ) {
                        error("(%f, %f) coordinates are different from (%.f, %f) at line %llu", pix_trs[i]->double_field_at(icols_x[i]), pix_trs[i]->double_field_at(icols_y[i]), mol_x, mol_y, nlines);
                    }
                }
                for(int32_t j=0; j < out_max_k; ++j) {
                    hprintf(wh_tsv, "\t%s", pix_trs[i]->str_field_at(icols_factor_K_idxs[i][j]));
                }
                for(int32_t j=0; j < out_max_p; ++j) {
                    hprintf(wh_tsv, "\t%s", pix_trs[i]->str_field_at(icols_factor_P_idxs[i][j]));
                }
            }
            hprintf(wh_tsv, "\n");
        }
        ++nlines;
        if ( nlines % 1000000 == 0 ) {
            notice("Processed %llu lines", nlines);
        }
    }
    hts_close(wh_tsv);
    notice("Processed %llu lines", nlines);
    notice("Analysis finished");

    return 0;
}
