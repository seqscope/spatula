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
int32_t cmdJoinPixelTSV(int32_t argc, char **argv)
{
    std::string in_mol_tsv;   // TSV file containing individual molecules
    std::vector<std::string> pix_prefix_tsvs; // vector of "[prefix],[tsv-path]" pairs for pixel-level projections
    std::string out_prefix;   // Output Prefix

    std::string out_suffix_tsv(".tsv.gz");            // suffix for the output TSV file
    std::string out_suffix_hist(".dist.hist.tsv");    // suffix for the histogram of match distance
    std::string out_suffix_summary(".summary.tsv");   // suffix for the summary file

    // output column names for tsv files
    //std::string colname_block("BLOCK");
    std::string colname_X("X");
    std::string colname_Y("Y");
    std::string sort_axis("X"); // both input files must be sorted by the same axis
    double bin_um = 1.0;        // unit to group temporary output
    double max_dist_um = 0.5;   // maximum distance between pixel-level output and the molecules 
    int32_t precision_um = 3;   // Output precision below the decimal point
    std::string csv_colnames_include; // comma-separated column names to include in the output TSV file
    std::string csv_colnames_exclude; // comma-separated column names to exclude in the output TSV file

    int32_t out_max_k = 1;            // maximum num of pixel-level factors to include in the joined output 
    int32_t out_max_p = 1;            // maximum num of pixel-level probs to include in the joined output
    double mu_scale = 1.0;            // scale factor for the mu values
    //bool skip_unmatched = false;      // do not print transcripts without matched factor

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("mol-tsv", &in_mol_tsv, "TSV file containing individual molecules")
    LONG_MULTI_STRING_PARAM("pix-prefix-tsv", &pix_prefix_tsvs, "TSV file containing pixel-level factors")
    LONG_STRING_PARAM("out-prefix", &out_prefix, "Output prefix for the joined TSV files")

    LONG_PARAM_GROUP("Key Parameters", NULL)
    LONG_DOUBLE_PARAM("bin-um", &bin_um, "Bin size for grouping the pixel-level output for indexing")
    LONG_DOUBLE_PARAM("max-dist-um", &max_dist_um, "Maximum distance in um to consider a match")
    LONG_DOUBLE_PARAM("mu-scale", &mu_scale, "Scale factor for the mu values")

    LONG_PARAM_GROUP("Expected columns in input and output", NULL)
    LONG_STRING_PARAM("colname-x", &colname_X, "Column name for X-axis")
    LONG_STRING_PARAM("colname-y", &colname_Y, "Column name for Y-axis")
    LONG_STRING_PARAM("sort-axis", &sort_axis, "Column name used in sorting. Both files must be sorted in the same axis (default: X-axis)")
    LONG_STRING_PARAM("colnames-include", &csv_colnames_include, "Comma-separated column names to include in the output TSV file")
    LONG_STRING_PARAM("colnames-exclude", &csv_colnames_exclude, "Comma-separated column names to exclude in the output TSV file")
    LONG_INT_PARAM("out-max-k", &out_max_k, "Maximum number of pixel-level factors to include in the joined output. (Default : 1)")
    LONG_INT_PARAM("out-max-p", &out_max_p, "Maximum number of pixel-level posterior probabilities to include in the joined output. (Default : 1)")

    LONG_PARAM_GROUP("Output File suffixes", NULL)
    LONG_STRING_PARAM("out-suffix-tsv", &out_suffix_tsv, "Suffix for the output TSV file")
    LONG_STRING_PARAM("out-suffix-hist", &out_suffix_hist, "Suffix for the histogram of match distance")
    LONG_STRING_PARAM("out-suffix-summary", &out_suffix_summary, "Suffix for the summary file")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check the input files
    if ( in_mol_tsv.empty() || pix_prefix_tsvs.empty() || out_prefix.empty() )
        error("--in-mol-tsv, --pix-prefix-tsv, and --out-prefix must be specified");

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

    bool is_sorted_by_x = false;
    if ( sort_axis.compare(colname_X) == 0 ) {
        is_sorted_by_x = true;
    }
    else if ( sort_axis.compare(colname_Y) == 0 ) {
        is_sorted_by_x = false;
    }
    else {
        error("Sorted axis %s matches to neither axis - %s and %s", sort_axis.c_str(), colname_X.c_str(), colname_Y.c_str());
    }

    //notice("is_sorted_by_x = %d", is_sorted_by_x);

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

    //notice("mol_colnames_set.size() = %zu", mol_colnames_set.size());

    // read the header and meta lines of pixel-level data
    std::vector<std::vector<std::string> > m_pix_colnames(n_pix);
    std::vector<double> v_offsets_x(n_pix);
    std::vector<double> v_offsets_y(n_pix);
    std::vector<double> v_scales(n_pix);
    std::vector<int32_t> v_topks(n_pix);
    //std::vector<int32_t> v_idxs_block(n_pix);
    std::vector<int32_t> v_idxs_x(n_pix);
    std::vector<int32_t> v_idxs_y(n_pix);
    std::vector<int32_t> v_idxs_k1(n_pix);
    std::vector<int32_t> v_idxs_p1(n_pix);
    for(int32_t i=0; i < n_pix; ++i) {
        //notice("i = %d", i);

        double offset_x = DBL_MAX;
        double offset_y = DBL_MAX;
        double scale = DBL_MAX;
        int32_t topk = 0;
        //bool block_axis_checked = false;
        int32_t idx_pix_x = -1, idx_pix_y = -1, idx_pix_k1 = -1, idx_pix_p1 = -1;
        tsv_reader* p_pix_tr = pix_trs[i];
        while( p_pix_tr->read_line() ) {
            const char* line = p_pix_tr->str_field_at(0);
            //notice("line = %s", line);
            if ( line[0] == '#' ) {
                if ( line[1] == '#' ) { // meta line, starting with '##'
                    std::string meta_line(p_pix_tr->str_field_at(0));
                    if ( p_pix_tr->nfields > 1 ) {
                        for(int i=1; i < p_pix_tr->nfields; ++i) {
                            meta_line += p_pix_tr->str_field_at(i);
                        }
                        line = meta_line.c_str(); // substitute the meta line removing blanks
                    }

                    std::vector<std::string> meta_toks;
                    split(meta_toks, ";", line+2);
                    for(int32_t i=0; i < (int32_t)meta_toks.size(); ++i) {
                        std::vector<std::string> keyvals;
                        split(keyvals, "=", meta_toks[i]);
                        if ( keyvals.size() != 2 ) {
                            error("Cannot parse '%s' in the meta line in %s. size = %zu, line = %s", meta_toks[i].c_str(), pix_tsvs[i].c_str(), keyvals.size(), line);
                        }
                        // if ( keyvals[0].compare("BLOCK_AXIS") == 0 ) {
                        //     // make sure that the axis is consistent
                        //     if ( keyvals[1].compare(sort_axis) != 0 ) {
                        //         error("BLOCK_AXIS=%s in %s, but sort_axis=%s", keyvals[1].c_str(), pix_tsvs[i].c_str(), sort_axis.c_str());
                        //     }
                        //     block_axis_checked = true;
                        // }
                        else if ( keyvals[0].compare("OFFSET_X") == 0 ) {
                            offset_x = atof(keyvals[1].c_str());
                        }
                        else if ( keyvals[0].compare("OFFSET_Y") == 0 ) {
                            offset_y = atof(keyvals[1].c_str());
                        }
                        else if ( keyvals[0].compare("SCALE") == 0 ) {
                            scale = atof(keyvals[1].c_str());
                        }
                        else if ( keyvals[0].compare("TOPK") == 0 ) {
                            topk = atoi(keyvals[1].c_str());
                        }
                    }
                }
                else { // header line, containing the column names
                    //notice("nfields = %d", p_pix_tr->nfields);
                    for(int32_t j=0; j < p_pix_tr->nfields; ++j) {
                        //notice("j=%d,  m_pix_colnames.size() = %zu, m_pix_colnames[i].size() = %zu", j, m_pix_colnames.size(), m_pix_colnames[i].size());
                        if ( j == 0 ) {
                            m_pix_colnames[i].push_back(line+1);
                        }
                        else {
                            m_pix_colnames[i].push_back(p_pix_tr->str_field_at(j));
                        }

                        std::string& colname = m_pix_colnames[i].back();
                        // if ( colname.compare(colname_block) == 0 ) {
                        //     idx_pix_block = j;
                        // }
                        if ( colname.compare(colname_X) == 0 ) {
                            idx_pix_x = j;
                        }
                        else if ( colname.compare(colname_Y) == 0 ) {
                            idx_pix_y = j;
                        }
                        else if ( colname.compare("K1") == 0 ) {
                            idx_pix_k1 = j;
                        }
                        else if ( colname.compare("P1") == 0 ) {
                            idx_pix_p1 = j;
                        }
                    }
                    break; // break if header line finishes parsing
                }
            }
            else {
                error("Non-header line found before the header line in the input TSV file %s", pix_tsvs[i].c_str());
            }
        }
        // sanity check to see whether required values are parsed
        if ( topk == 0 ) error("Cannot find TOPK field in the meta line of %s", pix_tsvs[i].c_str());
        if ( topk < out_max_k ) error("TOPK=%d is smaller than out-max-k=%d in %s", topk, out_max_k, pix_tsvs[i].c_str());
        if ( topk < out_max_p ) error("TOPK=%d is smaller than out-max-p=%d in %s", topk, out_max_p, pix_tsvs[i].c_str());
        if ( scale == DBL_MAX ) error("Cannot find SCALE field in the meta line of %s", pix_tsvs[i].c_str());
        if ( offset_x == DBL_MAX ) error("Cannot find OFFSET_X field in the meta line of %s", pix_tsvs[i].c_str());
        if ( offset_y == DBL_MAX ) error("Cannot find OFFSET_Y field in the meta line of %s", pix_tsvs[i].c_str());
        //if ( !block_axis_checked ) error("Cannot find BLOCK_AXIS field in the meta line of %s", pix_tsvs[i].c_str());
        //if ( idx_pix_block < 0 ) error("Cannot find %s in the header line of %s", colname_block.c_str(), pix_tsvs[i].c_str());
        if ( idx_pix_x < 0 ) error("Cannot find %s in the header line of %s", colname_X.c_str(), pix_tsvs[i].c_str());
        if ( idx_pix_y < 0 ) error("Cannot find %s in the header line of %s", colname_Y.c_str(), pix_tsvs[i].c_str());
        if ( idx_pix_k1 < 0 ) error("Cannot find K1 in the header line of %s", pix_tsvs[i].c_str());
        if ( idx_pix_p1 < 0 ) error("Cannot find P1 in the header line of %s", pix_tsvs[i].c_str());


        v_offsets_x[i] = offset_x;
        v_offsets_y[i] = offset_y;
        v_scales[i] = scale;
        v_topks[i] = topk;
        //v_idxs_block[i] = idx_pix_block;
        v_idxs_x[i] = idx_pix_x;
        v_idxs_y[i] = idx_pix_y;
        v_idxs_k1[i] = idx_pix_k1;
        v_idxs_p1[i] = idx_pix_p1;
    }


    //notice("foo1..");
    // read the header line of molecular-level data
    tsv_reader mol_tr(in_mol_tsv.c_str());
    std::vector<std::string> v_mol_colnames;
    int32_t idx_mol_x = -1, idx_mol_y = -1;
    if ( !mol_tr.read_line() ) {
        error("Cannot read the header line of %s", in_mol_tsv.c_str());
    }
    for(int32_t i=0; i < mol_tr.nfields; ++i) {
        v_mol_colnames.push_back(mol_tr.str_field_at(i));
        if ( colname_X.compare(v_mol_colnames.back()) == 0 ) {
            idx_mol_x = i;
        }
        else if ( colname_Y.compare(v_mol_colnames.back()) == 0 ) {
            idx_mol_y = i;
        }
    }
    if ( ( idx_mol_x < 0 ) || ( idx_mol_y < 0 ) ) {
        error("Cannot find %s or %s from the leader line of %s", colname_X.c_str(), colname_Y.c_str(), in_mol_tsv.c_str());
    }

    //notice("foo2..");

    // determine which columns to print (or not)
    std::vector<int32_t> v_mol_icols;
    for(int32_t i=0; i < v_mol_colnames.size(); ++i) {
        if ( mol_colnames_set.find(v_mol_colnames[i]) != mol_colnames_set.end() ) { // found the column name from the set
            if ( !mol_colnames_exclude ) {
                v_mol_icols.push_back(i);
            }
        }
        else if ( mol_colnames_exclude )  {
            v_mol_icols.push_back(i);
        }
    }
    if ( mol_colnames_exclude ) {
        if ( v_mol_icols.size() + mol_colnames_set.size() != v_mol_colnames.size() ) {
            error("Unrecognized or duplicate column names exist in --csv-mol-exclude parameter or %s: %s", csv_colnames_exclude.c_str(), in_mol_tsv.c_str());
        }
    }
    else {
        if ( v_mol_icols.size() != mol_colnames_set.size() ) {
            error("Unrecognized or duplicate columns in --csv-mol-include parameter or %s: %s", csv_colnames_include.c_str(), in_mol_tsv.c_str());
        }
    }

    //notice("foo3..");

    // construct simple bins to cache pixel-level factors
    std::vector<std::map<int32_t, std::map<int32_t, std::vector<pix_factor_t*> > > > v_bin2factors(n_pix);

    // read the first line of molecular-level data
    double mol_x, mol_y, mu_x, mu_y;
    std::vector<double> v_max_pix_vals(n_pix, -DBL_MAX);
    double max_major_val = -DBL_MAX;
    std::vector<std::map<int32_t, std::map<int32_t, std::vector<pix_factor_t*> > >::iterator > v_it(n_pix);
    std::vector<std::map<int32_t, std::vector<pix_factor_t*> >::iterator > v_it2(n_pix);

    // maintain the histogram of distances
    std::vector<std::map<int32_t, uint64_t> > v_dist2cnt(n_pix); // key: (int32_t)floor(log10(dist)*10), val : count

    //notice("Opening the output file..");
    // create output file
    htsFile* wh_tsv = hts_open((out_prefix + out_suffix_tsv).c_str(), out_suffix_tsv.compare(out_suffix_tsv.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");

    // write output header for the original TSV files
    for(int32_t i=0; i < v_mol_icols.size(); ++i) {
        if ( i > 0 ) {
            hprintf(wh_tsv, "\t");
        }
        hprintf(wh_tsv, "%s", v_mol_colnames[v_mol_icols[i]].c_str());
    }
    // hprintf(wh_tsv, "\tdist\tdiff\tncomp"); // DEBUG

    // write output header for each pixel-level file 
    for(int32_t i=0; i < n_pix; ++i) {
        for(int32_t j=0; j < out_max_k; ++j) {
            hprintf(wh_tsv, "\t%sK%d", pix_prefixes[i].c_str(), j+1);
        }
        for(int32_t j=0; j < out_max_p; ++j) {
            hprintf(wh_tsv, "\t%sP%d", pix_prefixes[i].c_str(), j+1);
        }
    }
    hprintf(wh_tsv, "\n");


    //notice("foo4..");

    notice("Started parsing the input file %s", in_mol_tsv.c_str());

    uint64_t n_full_match = 0, n_partial_match = 0, n_mol = 0, n_has_match = 0;
    std::vector<double> v_sums_match_dist(n_pix, 0.0);
    std::vector<int> v_n_matches(n_pix, 0); 
    
    while( mol_tr.read_line() ) {
        mol_x = mol_tr.double_field_at(idx_mol_x);
        mol_y = mol_tr.double_field_at(idx_mol_y);
        mu_x = mol_x / mu_scale; // micrometer-scale x value (consistent to the pixel-level data)
        mu_y = mol_y / mu_scale; // micrometer-scale y value (consistent to the pixel-level data)
        double new_max_major_val = is_sorted_by_x ? mu_x + max_dist_um : mu_y + max_dist_um;
        if ( new_max_major_val < max_major_val ) {
            error("Input file %s is not sorted by the axis %s", in_mol_tsv.c_str(), sort_axis.c_str());
        }
        max_major_val = new_max_major_val;

        // print the contents of TSV files
        for(int32_t i=0; i < v_mol_icols.size(); ++i) {
            if ( i > 0 ) {
                hprintf(wh_tsv, "\t");
            }
            hprintf(wh_tsv, "%s", mol_tr.str_field_at(v_mol_icols[i]));
        }

        //notice("mol_x = %lf, mol_y = %lf", mol_x, mol_y); 
        // read the pixel-level data up to the limit
        n_has_match = 0;
        for(int32_t i=0; i < n_pix; ++i) {
            uint64_t n_pix_read = 0;
            while ( ( v_max_pix_vals[i] < max_major_val ) && ( pix_trs[i]->read_line() ) ) { // keep reading the pixel-level data up to the limit
                // parse the line and fill in bin2factors
                pix_factor_t* ppf = new pix_factor_t( pix_trs[i], v_idxs_x[i], v_idxs_y[i], v_idxs_k1[i], v_idxs_p1[i], v_topks[i], v_offsets_x[i], v_offsets_y[i], v_scales[i]);
                int32_t bin_x = (int32_t)floor(ppf->x / bin_um);
                int32_t bin_y = (int32_t)floor(ppf->y / bin_um);
                int32_t bin_major = is_sorted_by_x ? bin_x : bin_y;
                int32_t bin_minor = is_sorted_by_x ? bin_y : bin_x;
                v_bin2factors[i][bin_major][bin_minor].push_back(ppf);
                if ( is_sorted_by_x ) {
                    if ( v_max_pix_vals[i] > ppf->x ) {
                        error("Pixel-level data in %s is not sorted by the axis %s", pix_tsvs[i].c_str(), sort_axis.c_str());
                    }
                    v_max_pix_vals[i] = ppf->x;
                }
                else {
                    if ( v_max_pix_vals[i] > ppf->y ) {
                        error("Pixel-level data in %s is not sorted by the axis %s", pix_tsvs[i].c_str(), sort_axis.c_str());
                    }
                    v_max_pix_vals[i] = ppf->y;
                }
                ++n_pix_read;
            }
            //if ( n_pix_read > 0 )
            //    notice("n_pix_read = %llu, max_major_val = %lf, max_pix_val = %lf", n_pix_read, max_major_val, max_pix_val);

            double min_major_val = is_sorted_by_x ? mu_x - max_dist_um : mu_y - max_dist_um;
            double min_minor_val = is_sorted_by_x ? mu_y - max_dist_um : mu_x - max_dist_um;
            double max_minor_val = is_sorted_by_x ? mu_y + max_dist_um : mu_x + max_dist_um;    

            int32_t min_major_bin = (int32_t)floor(min_major_val / bin_um);
            int32_t max_major_bin = (int32_t)floor(max_major_val / bin_um);

            // remove pixel-level data out of reach
            v_it[i] = v_bin2factors[i].begin();
            while( ( v_it[i]->first < min_major_bin ) && ( v_it[i] != v_bin2factors[i].end() ) ) {
                // remove all the pix_factor_t objects
                for(v_it2[i] = v_it[i]->second.begin(); v_it2[i] != v_it[i]->second.end(); ++v_it2[i]) {
                    for(int32_t j=0; j < (int32_t)v_it2[i]->second.size(); ++j) {
                        delete v_it2[i]->second[j];
                    }
                }
                ++v_it[i];
            }
            if ( v_it[i] == v_bin2factors[i].end() ) {
                v_bin2factors[i].clear();
            }
            else {
                v_bin2factors[i].erase(v_bin2factors[i].begin(), v_it[i]);
            }

            int32_t min_minor_bin = (int32_t)floor(min_minor_val / bin_um);
            int32_t max_minor_bin = (int32_t)floor(max_minor_val / bin_um);

            // identify best-matching bins to examine
            pix_factor_t* best_ppf = NULL;
            double best_dist = DBL_MAX; // max_dist_um * max_dist_um; 
            double next_dist = DBL_MAX; 
            double thres_dist = max_dist_um * max_dist_um;
            int32_t ncomp = 0;
            for(int32_t b1 = min_major_bin; b1 <= max_major_bin; ++b1) { // examine bins with limited range
                v_it[i] = v_bin2factors[i].find(b1);
                if ( v_it[i] != v_bin2factors[i].end() ) { // b1 exists
                    for(int32_t b2 = min_minor_bin; b2 <= max_minor_bin; ++b2) { // examine all minor axis
                        v_it2[i] = v_it[i]->second.find(b2);
                        if ( v_it2[i] != v_it[i]->second.end() ) { // b2 exists
                            std::vector<pix_factor_t*>& v = v_it2[i]->second;
                            // enumerate all pixel-level factors
                            for(int32_t j=0; j < (int32_t)v.size(); ++j) {
                                pix_factor_t* ppf = v[j];
                                double dist = (mu_x - ppf->x)*(mu_x - ppf->x) + (mu_y - ppf->y)*(mu_y - ppf->y);
                                ++ncomp;
                                if ( dist < best_dist ) {
                                    best_ppf = ppf;
                                    next_dist = best_dist;
                                    best_dist = dist;
                                }
                                else if ( dist < next_dist ) {
                                    next_dist = dist;
                                }
                            }
                        }
                    }
                }
            }

            if ( best_dist <= thres_dist ) {
                ++n_has_match;
                v_sums_match_dist[i] += sqrt(best_dist);
                ++v_n_matches[i];
            }

            for(int32_t j=0; j < out_max_k; ++j) {
                if ( best_dist > thres_dist ) {
                    hprintf(wh_tsv, "\tNA");
                }
                else {
                    hprintf(wh_tsv, "\t%u", (uint32_t)best_ppf->factors[j]);
                }
            }
            for(int32_t j=0; j < out_max_p; ++j) {
                if ( best_dist > thres_dist ) {
                    hprintf(wh_tsv, "\tNA");
                }
                else {
                    hprintf(wh_tsv, "\t%.3lg", best_ppf->probs[j]);
                }
            }
            //hprintf(wh_tsv, "\n");

            // update the histogram
            int32_t hist_key = (int32_t)floor(std::log10(best_dist)*10.0);
            ++v_dist2cnt[i][hist_key];
        }
        hprintf(wh_tsv, "\n");

        ++n_mol; // number of molecules processed
        if ( n_has_match > 0 ) {
            if ( n_has_match == n_pix ) {
                ++n_full_match;
            }
            else {
                ++n_partial_match;
            }
        }

        if ( n_mol % 1000000 == 0 ) {
            notice("Processed %llu transcripts and identified %llu full matches and %llu partial matches", n_mol, n_full_match, n_partial_match);
        }
    }
    hts_close(wh_tsv);
    notice("Processed %llu transcripts and identified %llu full matches and %llu partial matches", n_mol, n_full_match, n_partial_match);
    notice("Analysis finished");

    // Write the output histogram
    htsFile* wh_hist = hts_open((out_prefix + out_suffix_hist).c_str(), out_suffix_hist.compare(out_suffix_hist.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    hprintf(wh_hist, "index\tprefix\tmin\tmax\tcount\n");
    for(int32_t i=0; i < n_pix; ++i) {
        for(std::map<int32_t, uint64_t>::iterator it3 = v_dist2cnt[i].begin(); it3 != v_dist2cnt[i].end(); ++it3) {
            double min = pow(10, it3->first/10.0);
            double max = pow(10, (it3->first+1)/10.0);
            hprintf(wh_hist,"%d\t%s\t%.3lg\t%.3lg\t%llu\n", i, pix_prefixes[i].c_str(), min, max, it3->second);
        }
    }
    hts_close(wh_hist);

    return 0;
}
