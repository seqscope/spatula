#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

/////////////////////////////////////////////////////////////////////////
// convert-sge : Convert SGE format to a TSV format
////////////////////////////////////////////////////////////////////////
int32_t cmdSGE2TSV(int32_t argc, char **argv)
{
    std::string manifestf;
    std::string sgedir;
    std::string bcdf("barcodes.tsv.gz");
    std::string ftrf("features.tsv.gz");
    std::string mtxf("matrix.mtx.gz");
    std::string colname_gene_name("gene");
    std::string colname_gene_id("MoleculeID");
    std::string colname_count("Count");
    std::string colname_x("X");
    std::string colname_y("Y");
    bool print_feature_id = false;
    double units_per_um = 1.0;  // Output conversion factor 
    int32_t precision_um = 3;   // Output precision below the decimal point
    std::string outprefix;
    std::string out_suffix_tsv(".transcripts.tsv.gz");
    std::string out_suffix_ftr(".features.clean.tsv.gz");
    std::string out_suffix_minmax(".mimmax.tsv");
    int32_t icol_mtx = 1;
    std::string include_ftr_list;
    std::string exclude_ftr_list;
    std::string include_ftr_regex;
    std::string exclude_ftr_regex;
    std::string include_ftr_substr;
    std::string exclude_ftr_substr;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input Options", NULL)
    LONG_STRING_PARAM("sge", &sgedir, "SGE directory")
    LONG_STRING_PARAM("bcd", &bcdf, "Barcode file name")
    LONG_STRING_PARAM("ftr", &ftrf, "Feature file name")
    LONG_STRING_PARAM("mtx", &mtxf, "Matrix file name")
    LONG_STRING_PARAM("icol-mtx", &icol_mtx, "1-based column index in the matrix file to use as the count")

    LONG_PARAM_GROUP("Input Filtering Options", NULL)
    LONG_STRING_PARAM("include-feature-list", &include_ftr_list, "A file containing a list of input genes to be included (feature name of IDs)")
    LONG_STRING_PARAM("exclude-feature-list", &exclude_ftr_list, "A file containing a list of input genes to be excluded (feature name of IDs)")
    LONG_STRING_PARAM("include-feature-substr", &include_ftr_substr, "A substring of feature/gene names to be included")
    LONG_STRING_PARAM("exclude-feature-substr", &exclude_ftr_substr, "A substring of feature/gene names to be excluded")
    LONG_STRING_PARAM("include-feature-regex", &include_ftr_regex, "A regex pattern of feature/gene names to be included")
    LONG_STRING_PARAM("exclude-feature-regex", &exclude_ftr_regex, "A regex pattern of feature/gene names to be excluded")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_STRING_PARAM("out", &outprefix, "Output prefix")
    LONG_DOUBLE_PARAM("units-per-um", &units_per_um, "Coordinate unit per um")
    LONG_INT_PARAM("precision-um", &precision_um, "Output precision below the decimal point")

    LONG_PARAM_GROUP("Auxilary Output Options", NULL)
    LONG_PARAM("print-feature-id", &print_feature_id, "Print feature ID in output file")
    LONG_STRING_PARAM("colname-feature-name", &colname_gene_name, "Column name for feature/gene name")
    LONG_STRING_PARAM("colname-feature-id", &colname_gene_id, "Column name for feature/gene ID")
    LONG_STRING_PARAM("colname-count", &colname_count, "Column name for Count")
    LONG_STRING_PARAM("colname-x", &colname_x, "Column name for X")
    LONG_STRING_PARAM("colname-y", &colname_y, "Column name for Y")
    LONG_STRING_PARAM("suffix-tsv", &out_suffix_tsv, "Suffix for the transcript output file")
    LONG_STRING_PARAM("suffix-ftr", &out_suffix_ftr, "Suffix for the feature output file")
    LONG_STRING_PARAM("suffix-minmax", &out_suffix_minmax, "Suffix for the minmax output file")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (outprefix.empty() || sgedir.empty())
        error("--sge and --out must be specified");

    notice("Analysis started");

    // only one of the include and exclude options can be used
    int32_t include_sum = ( include_ftr_list.empty() ? 0 : 1 ) + ( include_ftr_regex.empty() ? 0 : 1 ) + ( include_ftr_substr.empty() ? 0 : 1 );
    if ( include_sum > 1 )
        error("Only one of --include-feature-list, --include-feature-regex, and --include-feature-substr can be used");
    
    int32_t exclude_sum = ( exclude_ftr_list.empty() ? 0 : 1 ) + ( exclude_ftr_regex.empty() ? 0 : 1 ) + ( exclude_ftr_substr.empty() ? 0 : 1 );
    if ( exclude_sum > 1 )
        error("Only one of --exclude-feature-list, --exclude-feature-regex, and --exclude-feature-substr can be used");

    if ( include_sum == 0 && exclude_sum == 0 )
        notice("No feature filtering is applied");
    else if ( include_sum > 0 && exclude_sum > 0 )
        warning("Both --include.. and --exclude.. filters applied. If both filters are matched, the feature will be excluded");

    // load the exclude and include gene lists
    std::set<std::string> include_ftr_set;
    std::set<std::string> exclude_ftr_set;
    if (!include_ftr_list.empty())
    {
        tsv_reader ftr_tr(include_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
            include_ftr_set.insert(ftr_tr.str_field_at(0));
    }
    if (!exclude_ftr_list.empty())
    {
        tsv_reader ftr_tr(exclude_ftr_list.c_str());
        while (ftr_tr.read_line() > 0)
            exclude_ftr_set.insert(ftr_tr.str_field_at(0));
    }

    // parse the gene list information
    tsv_reader ftr_tr((sgedir + "/" + ftrf).c_str());
    std::vector<std::string> ftr_ids;
    std::vector<std::string> ftr_names;
    std::vector<bool> ftr_passes;
    std::vector<uint64_t> ftr_cnts;
    std::regex regex_include(include_ftr_regex);
    std::regex regex_exclude(exclude_ftr_regex);
    while (ftr_tr.read_line() > 0)
    {
        // check if the gene is in the inclusion list
        bool include = true;
        if (!include_ftr_set.empty()) {
            if ( include_ftr_set.find(ftr_tr.str_field_at(0)) != include_ftr_set.end() ) {
                include = true;
            }
            else if ( include_ftr_set.find(ftr_tr.str_field_at(1)) != include_ftr_set.end() ) {
                include = true;
            }
            else {
                include = false;
            }
        }
        else if (!include_ftr_regex.empty()) {
            include = std::regex_search(ftr_tr.str_field_at(1), regex_include);
        }
        else if (!include_ftr_substr.empty()) {
            include = strstr(ftr_tr.str_field_at(1), include_ftr_substr.c_str()) != NULL;
        }

        bool exclude = false;
        if (!exclude_ftr_set.empty()) {
            if ( exclude_ftr_set.find(ftr_tr.str_field_at(0)) != exclude_ftr_set.end() ) {
                exclude = true;
            }
            else if ( exclude_ftr_set.find(ftr_tr.str_field_at(1)) != exclude_ftr_set.end() ) {
                exclude = true;
            }
            else {
                exclude = false;
            }
        }
        else if (!exclude_ftr_regex.empty()) {
            exclude = std::regex_search(ftr_tr.str_field_at(1), regex_exclude);
        }
        else if (!exclude_ftr_substr.empty()) {
            exclude = strstr(ftr_tr.str_field_at(1), exclude_ftr_substr.c_str()) != NULL;
        }

        ftr_ids.push_back(ftr_tr.str_field_at(0));
        ftr_names.push_back(ftr_tr.str_field_at(1));
        ftr_passes.push_back(include && !exclude);
        ftr_cnts.push_back(0);
    }

    // read the SGE matrix
    notice("Processing SGE directory %s ...", sgedir.c_str());
    sge_stream_reader ssr((sgedir + "/" + bcdf).c_str(), (sgedir + "/" + ftrf).c_str(), (sgedir + "/" + mtxf).c_str());
    htsFile* wh_tsv = hts_open((outprefix + out_suffix_tsv).c_str(), "wz");
    if ( wh_tsv == NULL )
        error("Failed to open %s%s for writing", outprefix.c_str(), out_suffix_tsv.c_str());
    // print the header
    if (print_feature_id)
        hprintf(wh_tsv, "%s\t%s\t%s\t%s\t%s\t%s\n", colname_x.c_str(), colname_y.c_str(),
            colname_gene_name.c_str(), colname_gene_id.c_str(), colname_count.c_str());
    else
        hprintf(wh_tsv, "%s\t%s\t%s\t%s\n", colname_x.c_str(), colname_y.c_str(),
            colname_gene_name.c_str(), colname_count.c_str());

    // read each line of the matrix file
    double um_x = 0.0, um_y = 0.0;
    double xmax = -1e99, ymax = -1e99, xmin = 1e99, ymin = 1e99;
    uint64_t intervals = 1000000;
    uint64_t nlines = 0;
    while (ssr.read_mtx())
    {
        if (nlines % intervals == 0)
            notice("Processing %llu lines... from %s/%s", nlines, sgedir.c_str(), mtxf.c_str());
        // remove lane and tile info and create global coordinates
        if (ssr.is_bcd_new)
        { // new barcode
            // calculate the new pixel coordinate
            um_x = ssr.cur_sbcd.px / units_per_um;
            um_y = ssr.cur_sbcd.py / units_per_um;

            // update the bounding box
            xmax = xmax > um_x ? xmax : um_x;
            ymax = ymax > um_y ? ymax : um_y;
            xmin = xmin < um_x ? xmin : um_x;
            ymin = ymin < um_y ? ymin : um_y;
        }
        // print individual transcripts
        if (ftr_passes[ssr.cur_iftr-1]) { // print only if the feature passes the filter
            //if ( floor(um_x) != um_x || floor(um_y) != um_y ) { // non-integer coordinates, use the precision
            hprintf(wh_tsv, "%.*f\t%.*f\t%s", precision_um, um_x, precision_um, um_y, ftr_names[ssr.cur_iftr-1].c_str());
            //}
            //else { // integer coordinates, ignore the precision
            //    hprintf(wh, "%d\t%d\t%s", int(um_x), int(um_y), ftr_names[ssr.cur_iftr-1].c_str());
            //}
            if (print_feature_id)
                hprintf(wh_tsv, "\t%s", ftr_ids[ssr.cur_iftr-1].c_str());

            // print the count information
            hprintf(wh_tsv, "\t%llu\n", ssr.cur_cnts[icol_mtx-1]);
        }
        ftr_cnts[ssr.cur_iftr-1] += ssr.cur_cnts[icol_mtx-1];
        ++nlines;
    }
    hts_close(wh_tsv);

    htsFile* wh_ftr = hts_open((outprefix + out_suffix_ftr).c_str(), "wz");
    if ( wh_ftr == NULL )
        error("Failed to open %s%s for writing", outprefix.c_str(), out_suffix_ftr.c_str());
    hprintf(wh_ftr, "%s\t%s\t%s\n", colname_gene_name.c_str(), colname_gene_id.c_str(), colname_count.c_str());
    for (int32_t i = 0; i < ftr_ids.size(); ++i)
    {
        if (ftr_passes[i] && ftr_cnts[i] > 0)
            hprintf(wh_ftr, "%s\t%s\t%llu\n", ftr_names[i].c_str(), ftr_ids[i].c_str(), ftr_cnts[i]);
    }
    hts_close(wh_ftr);

    htsFile* wh_minmax = hts_open((outprefix + out_suffix_minmax).c_str(), "w");
    if ( wh_minmax == NULL )
        error("Failed to open %s%s for writing", outprefix.c_str(), out_suffix_minmax.c_str());
    hprintf(wh_minmax, "xmin\t%.*f\n", precision_um, xmin);
    hprintf(wh_minmax, "xmax\t%.*f\n", precision_um, xmax);
    hprintf(wh_minmax, "ymin\t%.*f\n", precision_um, ymin);
    hprintf(wh_minmax, "ymax\t%.*f\n", precision_um, ymax);
    hts_close(wh_minmax);

    ssr.close();

    notice("Analysis finished");

    return 0;
}
