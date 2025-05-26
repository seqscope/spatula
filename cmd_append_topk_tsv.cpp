#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge2.h"
#include "file_utils.h"
#include "nlohmann/json.hpp"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <fstream>

/////////////////////////////////////////////////////////////////////////////////////////
// append-topk-tsv : Append TSV file by adding topK and topP columns
/////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdAppendTopKTSV(int32_t argc, char **argv)
{
    std::string in_tsv;   // TSV file containing individual molecules
    std::string in_json;  // JSON file containing the configuration of LDA4hex input
    std::string in_model; // Input model file (required with --reorder)
    std::string topK("topK"); // column name for topK
    std::string topP("topP"); // column name for topP
    std::string out_tsv;     // Output TSV file name
    std::string out_model;   // Output model file name (required with --reorder)
    bool reorder = false;    // reorder the columns in the output file
    bool keep_random_key = false; // keep the random key in the output file
    int32_t offset_model = 1;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("in-tsv", &in_tsv, "Input pseudobulk TSV file")
    LONG_STRING_PARAM("in-json", &in_json, "Input JSON file")
    LONG_STRING_PARAM("in-model", &in_model, "Input model file (required with --reorder)")
    LONG_STRING_PARAM("out-tsv", &out_tsv, "Output TSV file")
    LONG_STRING_PARAM("out-model", &out_model, "Output model file (required with --reorder)")
    LONG_PARAM("reorder", &reorder, "Reorder the columns in the output file based on the total count")
    LONG_PARAM("keep-random-key", &keep_random_key, "Keep the random key in the output file")

    LONG_PARAM_GROUP("Expected columns in input and output", NULL)
    LONG_INT_PARAM("offset-model", &offset_model, "Column index for the beginning of the input model file")
    LONG_STRING_PARAM("topK", &topK, "Column name for topK")
    LONG_STRING_PARAM("topP", &topK, "Column name for topP")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check the input files
    if ( in_tsv.empty() || out_tsv.empty() || in_json.empty() )
        error("--in-tsv, --in-json, and --out-tsv must be specified");

    // read the input JSON file
    std::ifstream json_file(in_json);
    if (!json_file.is_open()) {
        error("Cannot open JSON file %s", in_json.c_str());
    }
    nlohmann::json json_data;
    json_file >> json_data;

    // extract the header information
    int32_t icol_random_key = json_data["random_key"].get<int32_t>();
    int32_t offset_data = json_data["offset_data"].get<int32_t>();

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
        notice("icols.size() = %zu, offset_data = %d, offset_model = %d", icols.size(), offset_data, offset_model);
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
        //notice("icols.size() = %zu, offset_tsv = %d", icols.size(), offset_tsv);

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

        // write the reordered output model file
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
                for(int32_t i=offset_data; i < tr.nfields; ++i) {
                    if ( colnames[i-offset_data].compare(tr.str_field_at(i)) != 0 ) {
                        error("Incompatible input files between model and tsv file at column %d (%s vs %s)", i-offset_data, colnames[i-offset_data].c_str(), tr.str_field_at(i));
                    }
                }
            }
            else {
                // if the reorder was not requested, use the original order of the columns
                for(int32_t i=offset_data; i < tr.nfields; ++i) {
                    icols.push_back(i-offset_data);
                    colnames.push_back(tr.str_field_at(i));
                }
            }
            bool is_first_column = true;
            for(int32_t i=0; i < offset_data; ++i) {
                if ( !keep_random_key && ( i == icol_random_key ) )
                    continue;
                if ( !is_first_column ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(i));
                is_first_column = false;
            }
            for(int32_t i=0; i < (int32_t)icols.size(); ++i) {
                if ( !is_first_column ) hprintf(wh, "\t");
                hprintf(wh, "%s", colnames[i].c_str());
                is_first_column = false;
            }
            hprintf(wh, "\t%s\t%s\n", topK.c_str(), topP.c_str());
        }
        else { // process data lines
            double imax = 0;
            double maxP = tr.double_field_at(icols[0]+offset_data);
            for(int32_t i=1; i < (int32_t)icols.size(); ++i) {
                double p = tr.double_field_at(icols[i]+offset_data);
                if ( p > maxP ) {
                    maxP = p;
                    imax = i;
                }
            }
            bool is_first_column = true;
            for(int32_t i=0; i < offset_data; ++i) {
                if ( !keep_random_key && ( i == icol_random_key ) )
                    continue;
                if ( !is_first_column ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(i));
                is_first_column = false;
            }
            for(int32_t i=0; i < (int32_t)icols.size(); ++i) {
                if ( !is_first_column ) hprintf(wh, "\t");
                hprintf(wh, "%s", tr.str_field_at(icols[i]+offset_data));
                is_first_column = false;
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
