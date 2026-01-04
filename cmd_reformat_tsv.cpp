#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "generic_utils.h"
#include "sge.h"
#include <ctime>
#include <cmath>
#include <set>

struct refmt_scale_t {
public:
    double scale_factor;
    int32_t precision; // -1 means no precision control
    char format_type;  // 'f', 'e', 'g' are only possible value

    refmt_scale_t(double s, int32_t p, char t) : scale_factor(s), precision(p), format_type(t) {}
    refmt_scale_t() : scale_factor(1.0), precision(-1), format_type('f') {}

    refmt_scale_t(double s, const char* fmt) {
        scale_factor = s;
        if ( fmt[0] == '.' ) ++fmt; // ignore leading dot
        int32_t fmt_len = strlen(fmt);
        if ( fmt[fmt_len-1] == 'f' ) {
            format_type = 'f';
        }
        else if ( fmt[fmt_len-1] == 'e' ) {
            format_type = 'e';
        }
        else if ( fmt[fmt_len-1] == 'g' ) {
            format_type = 'g';
        }
        else {
            error("Unknown format type %s", fmt);
        }
        if ( fmt_len > 1 ) {
            precision = atoi(fmt);
        }
        else {
            precision = -1;
        }
    }

    void hprintf_scaled(htsFile* wf, double v) {
        v *= scale_factor;
        if ( precision < 0 ) {
            switch(format_type) {
                case 'f':
                    hprintf(wf, "%f", v);
                    break;
                case 'e':
                    hprintf(wf, "%e", v);
                    break;
                case 'g':
                    hprintf(wf, "%g", v);
                    break;
                default:
                    error("Unknown format type %c", format_type);
            }
        }
        else {
            switch(format_type) {
                case 'f':
                    hprintf(wf, "%.*f", precision, v);
                    break;
                case 'e':
                    hprintf(wf, "%.*e", precision, v);
                    break;
                case 'g':
                    hprintf(wf, "%.*g", precision, v);
                    break;
                default:
                    error("Unknown format type %c", format_type);
            }
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// reformat : Reformat TSV/CSV files by selecting or reordering columns
////////////////////////////////////////////////////////////////////////////////
int32_t cmdReformatTsv(int32_t argc, char **argv)
{
    std::string infile;
    std::string out_prefix;
    std::string in_delim("\t");
    std::string out_delim("\t");
    std::string colnames;
    std::vector<std::string> v_scales;
    int32_t skip_lines = 0;
    bool include_rest = false;
    bool unquote = false;
    bool write_minmax = false;
    bool write_features = false;
    bool add_count = false;
    std::string colname_x("X");
    std::string colname_y("Y");
    std::string colname_feature("gene");
    std::string colname_count("count");
    std::string suffix_tsv(".tsv.gz");
    std::string suffix_minmax(".minmax.tsv");
    std::string suffix_ftr(".features.tsv.gz");

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &infile, "Input CSV/TSV file")
    LONG_STRING_PARAM("in-delim", &in_delim, "Input delimiter")
    LONG_STRING_PARAM("colnames", &colnames, "Comma-separated column names to include in the output. To rename columns, use 'original_name:new_name' format. To add constants, use ':name:value' format.")
    LONG_MULTI_STRING_PARAM("scale", &v_scales, "Scale the columns by multiplying with the given value. Use 'colname:scale:format' format. use '0f' for integer, or use '3f', '3e', '5g', etc. Use the original column name before renaming")
    LONG_PARAM("include-rest", &include_rest, "Include all columns not specified in --colnames at the end")
    LONG_PARAM("unquote", &unquote, "Unquote the string values in the input file")
    LONG_INT_PARAM("skip-lines", &skip_lines, "Number of lines to skip at the beginning of the input file")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out_prefix, "Output CSV/TSV file prefix")
    LONG_STRING_PARAM("out-delim", &out_delim, "Output delimiter")
    LONG_STRING_PARAM("suffix-tsv", &suffix_tsv, "Output suffix for TSV files")
    LONG_STRING_PARAM("suffix-minmax", &suffix_minmax, "Output suffix for minmax files")
    LONG_STRING_PARAM("suffix-features", &suffix_ftr, "Output suffix for features files")
    LONG_STRING_PARAM("colname-x", &colname_x, "Column name for x-axis")
    LONG_STRING_PARAM("colname-y", &colname_y, "Column name for y-axis")
    LONG_STRING_PARAM("colname-feature", &colname_feature, "Column name for feature")
    LONG_STRING_PARAM("colname-count", &colname_count, "Column name for count")
    LONG_PARAM("write-minmax", &write_minmax, "Write minmax file")
    LONG_PARAM("write-features", &write_features, "Write features file")
    LONG_PARAM("add-count", &add_count, "Add count column to the output file, assigning 1 to all rows")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( infile.empty() || out_prefix.empty() ) 
        error("--in and --out must be specified");

    notice("Analysis started");

    tsv_reader tf(infile.c_str());
    if ( in_delim.size() != 1 ) {
        error("Input delimiter must be a single character");
    }
    tf.delimiter = in_delim.front();

    // skip the first few lines
    for (int32_t i = 0; i < skip_lines; ++i) {
        if ( !tf.read_line() ) {
            error("Cannot read the first %d lines from %s", skip_lines, infile.c_str());
        }
    }
    
    // read the header info
    if ( !tf.read_line() )
        error("Cannot read the header from %s", infile.c_str());

    // create a dictionary of input columns
    std::map<std::string, int32_t> col2idx;
    for(int32_t i=0; i < tf.nfields; ++i) {
        std::string s(tf.str_field_at(i));
        if ( unquote ) {
            unquote_string(s);
        }
        //notice("Column %d: %s", i, s.c_str());
        col2idx[s] = i;
    }

    // parse colname 
    std::vector<std::string> v_colnames;
    //std::vector<std::string> in_colnames;
    std::vector<std::string> out_colnames;
    std::vector<int32_t> out_colidx;
    std::set<int32_t> out_colidx_set;
    std::map<std::string,int32_t> out_col2idx;
    std::map<std::string,std::string> out_col2const;
    split(v_colnames, ",", colnames.c_str());
    for(int32_t i=0; i < (int32_t)v_colnames.size(); ++i) {
        std::string& colname = v_colnames[i];
        std::string new_colname;
        // if ':' is contained, rename the column
        if ( colname.front() == ':' ) { // add constant
            std::vector<std::string> toks;
            split(toks, ":", colname.substr(1).c_str());
            if ( toks.size() != 2 ) {
                error("Cannot parse %s", colname.c_str());
            }
            new_colname = toks[0];
            out_col2idx[new_colname] = -1;
            out_colnames.push_back(new_colname);
            out_colidx.push_back(-1);
            out_col2const[new_colname] = toks[1];
            continue;
        }
        else if ( colname.find(':') != std::string::npos ) {
            std::vector<std::string> toks;
            split(toks, ":", colname.c_str());
            if ( toks.size() != 2 ) {
                error("Cannot parse %s", colname.c_str());
            }
            colname = toks[0];
            new_colname = toks[1];
        }
        else {
            new_colname = colname;
        }

        int32_t idx = find_idx_by_key(col2idx, colname.c_str(), true);
        //in_colnames.push_back(colname);
        out_col2idx[new_colname] = idx;
        out_colnames.push_back(new_colname);
        out_colidx.push_back(idx);
        if ( out_colidx_set.insert(idx).second == false ) {
            error("Duplicated column name %s in %s", colname.c_str(), colnames.c_str());
        }
    }

    int32_t idx_x = -1, idx_y = -1, idx_feature = -1, idx_count = -1;
    double xmin = DBL_MAX, xmax = -DBL_MAX, ymin = DBL_MAX, ymax = -DBL_MAX;
    if ( write_minmax ) {
        if ( out_col2idx.find(colname_x) == out_col2idx.end() ) {
            error("Cannot find column %s in the input file. Needed when writing minmax file", colname_x.c_str());
        }
        idx_x = out_col2idx[colname_x];
        if ( out_col2idx.find(colname_y) == out_col2idx.end() ) {
            error("Cannot find column %s in the input file. Needed when writing minmax file", colname_y.c_str());
        }
        idx_y = out_col2idx[colname_y];
    }
    if ( write_features ) {
        if ( out_col2idx.find(colname_feature) == out_col2idx.end() ) {
            error("Cannot find column %s in the input file. Needed when writing features file", colname_feature.c_str());
        }
        idx_feature = out_col2idx[colname_feature];
        if ( out_col2idx.find(colname_count) == out_col2idx.end() ) {
            error("Cannot find column %s in the input file. Needed when writing features file", colname_count.c_str());
        }
        idx_count = out_col2idx[colname_count];
    }

    std::map<std::int32_t, refmt_scale_t> idx2scale;
    for (size_t i = 0; i < v_scales.size(); ++i) {
        std::string& scale = v_scales[i];
        std::vector<std::string> toks;
        split(toks, ":", scale.c_str());
        if ( ( toks.size() < 2 ) || ( toks.size() > 3 ) ) {
            error("Cannot parse %s", scale.c_str());
        }
        int32_t idx = find_idx_by_key(col2idx, toks[0].c_str(), true);
        double s = atof(toks[1].c_str());
        std::string fmt = (toks.size() == 3) ? toks[2] : "f";
        idx2scale[idx] = refmt_scale_t(s, fmt.c_str());
    }

    if ( include_rest ) {
        for(int32_t i=0; i < tf.nfields; ++i) {
            if ( out_colidx_set.find(i) == out_colidx_set.end() ) {
                out_colnames.push_back(tf.str_field_at(i));
                out_colidx.push_back(i);
            }
        }
    }

    int32_t expected_nfields = tf.nfields;

    std::string out_tsv = out_prefix + suffix_tsv;
    htsFile* wf = hts_open(out_tsv.c_str(), out_tsv.compare(out_tsv.size()-3, 3, ".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s", out_tsv.c_str());
    }

    // write the output header
    for (size_t i = 0; i < out_colnames.size(); ++i) {
        if ( i > 0 ) {
            hprintf(wf, "%s", out_delim.c_str());
        }
        hprintf(wf, "%s", out_colnames[i].c_str());
    }
    hprintf(wf, "\n");

    uint64_t nlines = 0;
    std::map<std::string, int32_t> ftr2cnt;
    while ( tf.read_line() ) {
        // write the output line
        if ( write_minmax ) {
            double x = tf.double_field_at(idx_x);
            double y = tf.double_field_at(idx_y);
            if ( x < xmin ) xmin = x;
            if ( x > xmax ) xmax = x;
            if ( y < ymin ) ymin = y;
            if ( y > ymax ) ymax = y;
        }
        if ( write_features ) {
            std::string ftr(tf.str_field_at(idx_feature));
            if ( unquote ) {
                unquote_string(ftr);
            }
            if ( idx_count < 0 ) {
                ++ftr2cnt[ftr];
            }
            else {
                ftr2cnt[ftr] += tf.int_field_at(idx_count);
            }
        }
        for (size_t i = 0; i < out_colidx.size(); ++i) {
            if ( i > 0 ) {
                hprintf(wf, "%s", out_delim.c_str());
            }
            if ( out_colidx[i] < 0 ) { // constant
                hprintf(wf, "%s", out_col2const[out_colnames[i]].c_str());
            }
            // if idx is in the scale list, print as scaled
            else if ( v_scales.size() > 0 && ( idx2scale.find(out_colidx[i]) != idx2scale.end() ) ) {
                double v = tf.double_field_at(out_colidx[i]);
                refmt_scale_t& scale = idx2scale[out_colidx[i]];
                scale.hprintf_scaled(wf, v);
            }
            else {
                std::string s(tf.str_field_at(out_colidx[i]));
                if ( unquote ) {
                    unquote_string(s);
                }
                hprintf(wf, "%s", s.c_str());
            }
        }
        hprintf(wf, "\n");
        if ( nlines % 1000000 == 0 ) {
            notice("Processing %llu lines...", nlines);
        }
        nlines++;
    }

    notice("Finished processing %llu lines and wrote %s", nlines, out_tsv.c_str());
    hts_close(wf);

    if ( write_minmax ) {
        std::string out_minmax = out_prefix + suffix_minmax;
        wf = hts_open(out_minmax.c_str(), out_minmax.compare(out_minmax.size()-3, 3, ".gz") == 0 ? "wz" : "w");
        if ( wf == NULL ) {
            error("Cannot open output file %s", out_minmax.c_str());
        }
        hprintf(wf, "xmin\t%f\n", xmin);
        hprintf(wf, "xmax\t%f\n", xmax);
        hprintf(wf, "ymin\t%f\n", ymin);
        hprintf(wf, "ymax\t%f\n", ymax);
        hts_close(wf);
        notice("Finished writing minmax file %s", out_minmax.c_str());
    }

    if ( write_features ) {
        std::string out_ftr = out_prefix + suffix_ftr;
        wf = hts_open(out_ftr.c_str(), out_ftr.compare(out_ftr.size()-3, 3, ".gz") == 0 ? "wz" : "w");
        if ( wf == NULL ) {
            error("Cannot open output file %s", out_ftr.c_str());
        }
        //hprintf(wf, "%s\t%s\n", colname_feature.c_str(), colname_count.c_str());
        for (auto it = ftr2cnt.begin(); it != ftr2cnt.end(); ++it) {
            hprintf(wf, "%s\t%llu\n", it->first.c_str(), it->second);
        }
        hts_close(wf);
        notice("Finished writing features file %s", out_ftr.c_str());
    }

    notice("Analysis finished");


    return 0;
}
