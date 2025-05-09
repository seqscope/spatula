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
    std::string in_model; // Input model file (required with --reorder)
    std::string topK("topK"); // column name for topK
    std::string topP("topP"); // column name for topP
    std::string out_tsv;     // Output TSV file name
    std::string out_model;   // Output model file name (required with --reorder)
    bool reorder = false;    // reorder the columns in the output file
    int32_t offset_tsv = 2;
    int32_t offset_model = 1;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("in-tsv", &in_tsv, "Input pseudobulk TSV file")
    LONG_STRING_PARAM("in-model", &in_model, "Input model file (required with --reorder)")
    LONG_STRING_PARAM("out-tsv", &out_tsv, "Output TSV file")
    LONG_STRING_PARAM("out-model", &out_model, "Output model file (required with --reorder)")
    LONG_PARAM("reorder", &reorder, "Reorder the columns in the output file based on the total count")

    LONG_PARAM_GROUP("Expected columns in input and output", NULL)
    LONG_INT_PARAM("offset-model", &offset_model, "Column index for the beginning of the input model file")
    LONG_INT_PARAM("offset-tsv", &offset_tsv, "Column index for the beginning of the input TSV file")
    LONG_STRING_PARAM("topK", &topK, "Column name for topK")
    LONG_STRING_PARAM("topP", &topK, "Column name for topP")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check the input files
    if ( in_tsv.empty() || out_tsv.empty() )
        error("--in-tsv, and --out-tsv must be specified");

    std::vector<int32_t> icols;
    std::vector<std::string> colnames;
    uint64_t nlines = 0;
    if ( reorder ) {
        if ( in_model.empty() ) {
            error("--in-model must be specified with --reorder");
        }
        if ( out_model.empty() ) {
            error("--out-model must be specified with --reorder");
        }
        // read the model file
        tsv_reader tr(in_model.c_str());
        std::vector<double> v_colsums;
        notice("icols.size() = %zu, offset_tsv = %d", icols.size(), offset_tsv);
        while( tr.read_line() ) {
            if ( nlines == 0 ) {
                notice("tr.nfields = %d", tr.nfields);
                for(int32_t i=offset_model; i < tr.nfields; ++i) {
                    icols.push_back(i-offset_model);
                    colnames.push_back(tr.str_field_at(i));
                    v_colsums.push_back(0.0);
                }
            }
            else { // process data lines
                for(int32_t i=0; i < (int32_t)icols.size(); ++i) {
                    double p = tr.double_field_at(icols[i]+offset_model);
                    v_colsums[i] += p;
                }
            }
            ++nlines;
        }
        tr.close();
        notice("icols.size() = %zu, offset_tsv = %d", icols.size(), offset_tsv);

        // sort the columns based on the total count
        std::vector<std::pair<double, int32_t> > v_colsums_idx;
        for(int32_t i=0; i < (int32_t)icols.size(); ++i) {
            v_colsums_idx.push_back(std::make_pair(v_colsums[i], icols[i]));
        }
        std::sort(v_colsums_idx.begin(), v_colsums_idx.end(), std::greater<std::pair<double, int32_t> >());

        // reorder the columns
        icols.clear();
        for(int32_t i=0; i < (int32_t)v_colsums_idx.size(); ++i) {
            icols.push_back(v_colsums_idx[i].second);
            //notice("%d\t%d", i, icols[i]);
        }

        // write the output model file
        tsv_reader tr2(in_model.c_str());
        htsFile* wh = hts_open(out_model.c_str(), out_model.compare(out_model.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
        nlines = 0;
        while( tr2.read_line() ) {
            for(int32_t i=0; i < offset_model; ++i) {
                if ( i > 0 ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr2.str_field_at(i));
            }
            for(int32_t i=0; i < (int32_t)icols.size(); ++i) {
                if ( i + offset_model > 0 ) hprintf(wh, "\t");
                if ( nlines == 0 ) {
                    hprintf(wh, "%s", colnames[i].c_str());
                }
                else {
                    hprintf(wh, "%s", tr2.str_field_at(icols[i]+offset_model));
                }
            }
            hprintf(wh, "\n");
            ++nlines;
        }
        tr2.close();
        hts_close(wh);
        notice("Reordered the columns in the output model file %s", out_model.c_str());
    }
    else {
        notice("Reordering is not requested, so the output model file will not be generated");
    }

    tsv_reader tr(in_tsv.c_str());
    htsFile* wh = hts_open(out_tsv.c_str(), out_tsv.compare(out_tsv.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");

    nlines = 0;
    while( tr.read_line() ) {
        if ( nlines == 0 ) { // process header
            if ( reorder ) {
                for(int32_t i=offset_tsv; i < tr.nfields; ++i) {
                    if ( colnames[i-offset_tsv].compare(tr.str_field_at(i)) != 0 ) {
                        error("Incompatible input files between model and tsv file at column %d (%s vs %s)", i-offset_tsv, colnames[i-offset_tsv].c_str(), tr.str_field_at(i));
                    }
                }
            }
            else {
                for(int32_t i=offset_tsv; i < tr.nfields; ++i) {
                    icols.push_back(i-offset_tsv);
                    colnames.push_back(tr.str_field_at(i));
                }
            }
            for(int32_t i=0; i < offset_tsv; ++i) {
                if ( i > 0 ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(i));
            }
            for(int32_t i=0; i < (int32_t)icols.size(); ++i) {
                if ( i + offset_tsv > 0 ) hprintf(wh, "\t");
                //hprintf(wh, "%s", tr.str_field_at(icols[i]+offset_tsv));
                hprintf(wh, "%s", colnames[i].c_str());
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
            for(int32_t i=0; i < offset_tsv; ++i) {
                if ( i > 0 ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(i));
            }
            for(int32_t i=0; i < (int32_t)icols.size(); ++i) {
                if ( i + offset_tsv > 0 ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(icols[i]+offset_tsv));
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
