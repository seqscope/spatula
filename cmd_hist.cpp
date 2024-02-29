#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include <ctime>
#include <cmath>
#include <set>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

std::string double2decstr(double original, int32_t digits) { 
    char buf[255];
    snprintf(buf, 255, "%.*f", digits, original);
    return std::string(buf);
}
std::string double2unitstr(double original, double unit) {
    char buf[255];
    int32_t exponent_o = floor(log10(abs(original)));    
    if ( ( original == floor(original) ) && ( exponent_o < 5 ) ) {
        snprintf(buf, 255, "%d", (int32_t)original);
    }
    else {
        int32_t exponent_u = floor(log10(abs(unit)));
        int32_t sig_digits = exponent_o - exponent_u + 1;
        if ( sig_digits < 2 ) sig_digits = 2; 
        snprintf(buf, 255, "%.*g", sig_digits, original);
    }
    return std::string(buf);   
}

double double_approx(double original, int32_t sig_digits = 0) {
    if ( original == 0 ) return 0;
    int32_t exponent = floor(log10(abs(original)));
    if ( sig_digits == 0 ) sig_digits = exponent;

    double significand = original / pow(10.0, exponent);
    double factor = pow(10.0, sig_digits);
    significand = round(significand * factor) / factor;
    double approx = significand * pow(10.0, exponent);
    return approx;
}

void find_minmax(const std::vector<double>& vals, double& min_x, double& max_x) {
    min_x = vals[0];
    max_x = vals[0];
    for(int32_t i=1; i < (int32_t)vals.size(); ++i) {
        if ( vals[i] < min_x ) min_x = vals[i];
        if ( vals[i] > max_x ) max_x = vals[i];
    }
}

/////////////////////////////////////////////////////////////////////////
// hist : Draw a text-based histogram based on input data 
////////////////////////////////////////////////////////////////////////
int32_t cmdHist(int32_t argc, char **argv)
{
    std::string tsvf("-");
    std::string outf("-");
    int32_t icol = 1;
    int32_t batch_size = 1000000;
    double bin_width = 0.0;
    int32_t num_bins = 0;
    bool show_fraction = false;
    bool show_cumulative = false;
    bool show_median = false;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("file", &tsvf, "Input file. Use - for stdin")
    LONG_INT_PARAM("column", &icol, "1-based index of the column to select")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output tsv file to show the results. Use - for stdout")
    LONG_DOUBLE_PARAM("bin-width", &bin_width, "Width of the each bin")
    LONG_INT_PARAM("num-bins", &num_bins, "Number of bins (10 by default)")
    LONG_PARAM("show-fraction", &show_fraction, "Show the fraction of the total")
    LONG_PARAM("show-cumulative", &show_cumulative, "Show the cumulative fraction")
    LONG_PARAM("show-median", &show_median, "Show the median value for each interval")

    LONG_PARAM_GROUP("Other settings", NULL)
    LONG_INT_PARAM("batch-size", &batch_size, "Size of initial batch to determine the bin width")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( tsvf.empty() || outf.empty() )
        error("--tsv and --out must be specified");

    if ( num_bins > 0 && bin_width > 0 )
        error("--num-bins and --bin-width cannot be specified at the same time");
    
    if ( num_bins == 0 && bin_width == 0 )
        num_bins = 10;

    notice("Analysis started");

    tsv_reader tf(tsvf.c_str());
    std::vector<double> batch;
    std::map<int64_t, uint64_t> hist;
    uint64_t nlines = 0;
    double min_x, max_x;

    while ( tf.read_line() ) {
        if ( tf.nfields < icol ) 
            error("Input file %s does not have enough columns - only %d", tsvf.c_str(), tf.nfields);

        double x = tf.double_field_at(icol-1);
        if ( bin_width == 0 ) {
            batch.push_back(x);
            if ( batch.size() >= batch_size ) {
                find_minmax(batch, min_x, max_x);
                bin_width = double_approx((max_x - min_x) / num_bins, 1);
                notice("Determined the bin width = %lf", bin_width);
                for(int32_t i=0; i < (int32_t)batch.size(); ++i) {
                    int64_t bin = (int64_t)(batch[i] / bin_width);
                    ++hist[bin];
                }
                batch.clear();
            }
        }
        else {
            int64_t bin = (int64_t)(x / bin_width);
            ++hist[bin];
        }
        ++nlines;
    }

    if ( nlines == 0 )
        error("No data found in the input file");

    if ( bin_width == 0 ) {
        find_minmax(batch, min_x, max_x);
        bin_width = double_approx((max_x - min_x) / num_bins, 1);
        notice("Determined the bin width = %lf", bin_width);
        for(int32_t i=0; i < (int32_t)batch.size(); ++i) {
            int64_t bin = (int64_t)(batch[i] / bin_width);
            ++hist[bin];
        }
        batch.clear();
    }

    // print the histogram
    htsFile* wh = hts_open(outf.c_str(), "w");
    if ( wh == NULL )
        error("Cannot open the output file %s", outf.c_str());
    
    std::vector<std::string> header;
    header.push_back("from");
    header.push_back("to");
    if ( show_median ) header.push_back("median");
    header.push_back("count");
    if ( show_fraction ) header.push_back("frac");
    if ( show_cumulative ) {
        header.push_back("cumul");
        if ( show_fraction )
            header.push_back("fcumul");
    }
    hprintf(wh, "%s", header[0].c_str());
    for(int32_t i=1; i < (int32_t)header.size(); ++i)
        hprintf(wh, "\t%s", header[i].c_str());
    hprintf(wh, "\n");

    
    std::map<int64_t, uint64_t>::iterator it;
    uint64_t cumul = 0;
    for(it = hist.begin(); it != hist.end(); ++it) {
        double from = it->first * bin_width;
        double to = from + bin_width;
        std::vector<std::string> row;
        row.push_back(double2unitstr(from, bin_width));
        row.push_back(double2unitstr(to, bin_width));
        if ( show_median ) {
            double median = (from + to) / 2.0;
            row.push_back(double2unitstr(median, bin_width/10));
        }
        row.push_back(std::to_string(it->second));
        // print the fraction
        if ( show_fraction ) 
            row.push_back(double2decstr((double)it->second / nlines,5));
        // print the cumulative fraction
        if ( show_cumulative ) {
            cumul += it->second;
            row.push_back(std::to_string(cumul));
            if ( show_fraction )
                row.push_back(double2decstr((double)cumul / nlines,5));
        }
        hprintf(wh, "%s", row[0].c_str());
        for(int32_t i=1; i < (int32_t)row.size(); ++i)
            hprintf(wh, "\t%s", row[i].c_str());
        hprintf(wh, "\n");
    }
    hts_close(wh);

    notice("Analysis started");

    return 0;
}
