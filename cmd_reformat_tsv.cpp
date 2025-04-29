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
    std::string outfile;
    std::string in_delim("\t");
    std::string out_delim("\t");
    std::string colnames;
    std::vector<std::string> v_scales;
    int32_t skip_lines = 0;
    bool include_rest = false;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &infile, "Input CSV/TSV file")
    LONG_STRING_PARAM("in-delim", &in_delim, "Input delimiter")
    LONG_STRING_PARAM("colnames", &colnames, "Comma-separated column names to include in the output. To rename columns, use 'original_name:new_name' format")
    LONG_MULTI_STRING_PARAM("scale", &v_scales, "Scale the columns by multiplying with the given value. Use 'colname:scale:format' format. use '0f' for integer, or use '3f', '3e', '5g', etc. Use the original column name before renaming")
    LONG_PARAM("include-rest", &include_rest, "Include all columns not specified in --colnames at the end")
    LONG_PARAM("skip-lines", &skip_lines, "Number of lines to skip at the beginning of the input file")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outfile, "Output CSV/TSV file")
    LONG_STRING_PARAM("out-delim", &out_delim, "Output delimiter")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( infile.empty() || outfile.empty() ) 
        error("--in and --out must be specified");

    notice("Analysis started");

    tsv_reader tf(infile.c_str());

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
        const char* s = tf.str_field_at(i);
        col2idx[s] = i;
    }

    // parse colname 
    std::vector<std::string> v_colnames;
    std::vector<std::string> in_colnames;
    std::vector<std::string> out_colnames;
    std::vector<int32_t> out_colidx;
    std::set<int32_t> out_colidx_set;
    split(v_colnames, ",", colnames.c_str());
    for(int32_t i=0; i < (int32_t)v_colnames.size(); ++i) {
        std::string& colname = v_colnames[i];
        std::string new_colname;
        // if ':' is contained, rename the column
        if ( colname.find(':') != std::string::npos ) {
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
        in_colnames.push_back(colname);
        out_colnames.push_back(new_colname);
        out_colidx.push_back(idx);
        if ( out_colidx_set.insert(idx).second == false ) {
            error("Duplicated column name %s in %s", colname.c_str(), colnames.c_str());
        }
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

    htsFile* wf = hts_open(outfile.c_str(), outfile.compare(outfile.size()-3, 3, ".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s", outfile.c_str());
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
    while ( tf.read_line() ) {
        // write the output line
        for (size_t i = 0; i < out_colidx.size(); ++i) {
            if ( i > 0 ) {
                hprintf(wf, "%s", out_delim.c_str());
            }
            // if idx is in the scale list, print as scaled
            if ( v_scales.size() > 0 && ( idx2scale.find(out_colidx[i]) != idx2scale.end() ) ) {
                double v = tf.double_field_at(out_colidx[i]);
                refmt_scale_t& scale = idx2scale[out_colidx[i]];
                scale.hprintf_scaled(wf, v);
            }
            else {
                const char* s = tf.str_field_at(out_colidx[i]);
                hprintf(wf, "%s", s);
            }
        }
        hprintf(wf, "\n");
        if ( nlines % 1000000 == 0 ) {
            notice("Processing %llu lines...", nlines);
        }
        nlines++;
    }

    notice("Finished processing %llu lines", nlines);
    notice("Analysis finished");

    hts_close(wf);

    return 0;
}
