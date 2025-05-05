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

/////////////////////////////////////////////////////////////////////////////////////////
// append-topk-tsv : Append TSV file by adding topK and topP columns
/////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdAppendTopKTSV(int32_t argc, char **argv)
{
    std::string in_tsv;   // TSV file containing individual molecules
    std::string topK("topK"); // column name for topK
    std::string topP("topP"); // column name for topP
    std::string out_tsv;     // Output TSV file name
    int32_t icol_beg = 2;    // Column index for the beginning of the input TSV file
    int32_t icol_end = -1;   // Column index for the end of the input TSV file

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("in-tsv", &in_tsv, "Input TSV file")
    LONG_STRING_PARAM("out-tsv", &out_tsv, "Output TSV file")

    LONG_PARAM_GROUP("Expected columns in input and output", NULL)
    LONG_INT_PARAM("icol-beg", &icol_beg, "Column index for the beginning of the input TSV file")
    LONG_INT_PARAM("icol-end", &icol_end, "Column index for the end of the input TSV file")
    LONG_STRING_PARAM("topK", &topK, "Column name for topK")
    LONG_STRING_PARAM("topP", &topK, "Column name for topP")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check the input files
    if ( in_tsv.empty() || out_tsv.empty() )
        error("--in-tsv, and --out-tsv must be specified");

    tsv_reader tr(in_tsv.c_str());
    htsFile* wh = hts_open(out_tsv.c_str(), out_tsv.compare(out_tsv.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");

    uint64_t nlines = 0;
    std::vector<int32_t> icols;
    std::vector<std::string> colnames;
    while( tr.read_line() ) {
        if ( nlines == 0 ) { // process header
            if ( icol_end < 0 ) icol_end = tr.nfields - 1;
            for(int32_t i=icol_beg; i <= icol_end; ++i) {
                icols.push_back(i);
                colnames.push_back(tr.str_field_at(i));
            }
            for(int32_t i=0; i < tr.nfields; ++i) {
                if ( i > 0 ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(i));
            }
            hprintf(wh, "\t%s\t%s\n", topK.c_str(), topP.c_str());
        }
        else { // process data lines
            double imax = 0;
            double maxP = tr.double_field_at(icols[0]);
            for(int32_t i=1; i < (int32_t)icols.size(); ++i) {
                double p = tr.double_field_at(icols[i]);
                if ( p > maxP ) {
                    maxP = p;
                    imax = i;
                }
            }
            for(int32_t i=0; i < tr.nfields; ++i) {
                if ( i > 0 ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(i));
            }
            hprintf(wh, "\t%s\t%.5g\n", colnames[imax].c_str(), maxP);
        }
        ++nlines;
        if ( nlines % 1000000 == 0 ) {
            notice("Processed %llu lines", nlines);
        }
    }
    notice("Finished processing %llu lines", nlines);

    hts_close(wh);
    tr.close();
    notice("Analysis finished");

    return 0;
}
