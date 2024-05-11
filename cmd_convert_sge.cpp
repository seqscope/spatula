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

/////////////////////////////////////////////////////////////////////////
// convert-sge : Convert SGE format to a TSV format
////////////////////////////////////////////////////////////////////////
int32_t cmdConvertSGE(int32_t argc, char **argv)
{
    // core input files
    std::string in_sgedir;
    std::string outdir;
    bool out_sge = false;
    bool out_tsv = false;
    std::string posf;   // external position files

    // SGE and TSV file names
    std::string sgedir("sge");
    std::string tsvdir("tsv");
    std::string sge_bcdf("barcodes.tsv.gz");
    std::string sge_ftrf("features.tsv.gz");
    std::string sge_mtxf("matrix.mtx.gz");
    std::string tsv_mtxf("transcripts.tsv.gz");
    std::string tsv_ftrf("features.clean.tsv.gz");
    std::string tsv_minmaxf("minmax.tsv");

    // input column indices
    int32_t in_icol_bcd_barcode = 1;
    int32_t in_icol_bcd_px = 6; // 1-based column index in the barcode file to use as the X coordinate
    int32_t in_icol_bcd_py = 7; // 1-based column index in the barcode file to use as the Y coordinate
    int32_t in_icol_ftr_id = 1;
    int32_t in_icol_ftr_name = 2;
    int32_t in_icol_mtx = 1;

    // output column names for tsv files
    std::string colname_gene_name("gene");
    std::string colname_gene_id("MoleculeID");
    std::string colname_count("Count");
    std::string colname_x("X");
    std::string colname_y("Y");
    bool print_feature_id = false;

    // common output options
    double units_per_um = 1.0;  // Output conversion factor 
    int32_t precision_um = 3;   // Output precision below the decimal point
    std::string include_ftr_list;
    std::string exclude_ftr_list;
    std::string include_ftr_regex;
    std::string exclude_ftr_regex;
    std::string include_ftr_substr;
    std::string exclude_ftr_substr;

    // additional barcode position file
    std::string pos_colname_barcode("barcode");
    std::string pos_colname_x("pxl_row_in_fullres");
    std::string pos_colname_y("pxl_col_in_fullres");
    std::string pos_delim(",");

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("in-sge", &in_sgedir, "Input SGE directory")
    LONG_STRING_PARAM("out-dir", &outdir, "Output directory")
    LONG_PARAM("out-sge", &out_sge, "Generate SGE output")
    LONG_PARAM("out-tsv", &out_tsv, "Generate TSV output")

    LONG_PARAM_GROUP("File names for SGE/TSV input/output", NULL)
    LONG_STRING_PARAM("sge-bcd", &sge_bcdf, "Barcode file name in SGE directory")
    LONG_STRING_PARAM("sge-ftr", &sge_ftrf, "Feature file name in SGE directory")
    LONG_STRING_PARAM("sge-mtx", &sge_mtxf, "Matrix file name in SGE directory")
    LONG_STRING_PARAM("tsv-mtx", &tsv_mtxf, "Transcript file name in TSV directory")
    LONG_STRING_PARAM("tsv-ftr", &tsv_ftrf, "Feature file name in TSV directory")
    LONG_STRING_PARAM("tsv-minmax", &tsv_minmaxf, "Minmax file name in TSV directory")

    LONG_PARAM_GROUP("Expected column index in SGE input", NULL)
    LONG_INT_PARAM("icol-mtx", &in_icol_mtx, "1-based column index in the matrix file to use as the count")
    LONG_INT_PARAM("icol-bcd-barcode", &in_icol_bcd_barcode, "1-based column index of barcode in the barcode file")
    LONG_INT_PARAM("icol-bcd-x", &in_icol_bcd_px, "1-based column index of x coordinate in the barcode file")
    LONG_INT_PARAM("icol-bcd-y", &in_icol_bcd_py, "1-based column index of y coordinate in the barcode file")
    LONG_INT_PARAM("icol-ftr-id", &in_icol_ftr_id, "1-based column index of feature ID in the barcode file")
    LONG_INT_PARAM("icol-ftr-name", &in_icol_ftr_name, "1-based column index of feature name in the barcode file")

    LONG_PARAM_GROUP("Additional Barcode Position File", NULL)
    LONG_STRING_PARAM("pos", &posf, "Position file name that contains separate X and Y coordinates")
    LONG_STRING_PARAM("pos-colname-barcode", &pos_colname_barcode, "Column name for barcode in the position file")
    LONG_STRING_PARAM("pos-colname-x", &pos_colname_x, "Column name for X-axis in the position file")
    LONG_STRING_PARAM("pos-colname-y", &pos_colname_y, "Column name for Y-axis in the position file")
    LONG_STRING_PARAM("pos-delim", &pos_delim, "Delimiter for the position file (default: \",\")")

    LONG_PARAM_GROUP("Input Filtering Options", NULL)
    LONG_STRING_PARAM("include-feature-list", &include_ftr_list, "A file containing a list of input genes to be included (feature name of IDs)")
    LONG_STRING_PARAM("exclude-feature-list", &exclude_ftr_list, "A file containing a list of input genes to be excluded (feature name of IDs)")
    LONG_STRING_PARAM("include-feature-substr", &include_ftr_substr, "A substring of feature/gene names to be included")
    LONG_STRING_PARAM("exclude-feature-substr", &exclude_ftr_substr, "A substring of feature/gene names to be excluded")
    LONG_STRING_PARAM("include-feature-regex", &include_ftr_regex, "A regex pattern of feature/gene names to be included")
    LONG_STRING_PARAM("exclude-feature-regex", &exclude_ftr_regex, "A regex pattern of feature/gene names to be excluded")

    LONG_PARAM_GROUP("Key Output Options", NULL)
    LONG_DOUBLE_PARAM("units-per-um", &units_per_um, "Coordinate unit per um (conversion factor)")
    LONG_INT_PARAM("precision-um", &precision_um, "Output precision below the decimal point")

    LONG_PARAM_GROUP("Auxilary Output Options for TSV output", NULL)
    LONG_PARAM("print-feature-id", &print_feature_id, "Print feature ID in output file")
    LONG_STRING_PARAM("colname-feature-name", &colname_gene_name, "Column name for feature/gene name")
    LONG_STRING_PARAM("colname-feature-id", &colname_gene_id, "Column name for feature/gene ID")
    LONG_STRING_PARAM("colname-count", &colname_count, "Column name for Count")
    LONG_STRING_PARAM("colname-x", &colname_x, "Column name for X")
    LONG_STRING_PARAM("colname-y", &colname_y, "Column name for Y")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( in_sgedir.empty() && outdir.empty() )
        error("Input/output directory --in-sge and --out-dir must be specified");

    if ( !( out_sge || out_tsv ) )
        error("At least one of --out-sge or --out-tsv must specified");


    // remove trailing slash if exists
    if ( outdir[outdir.size()-1] == '/' )
        outdir = outdir.substr(0, outdir.size()-1);

    if ( out_sge ) { // create sge directory
        makePath((outdir + "/" + sgedir).c_str());
    }
    if ( out_tsv ) { // create tsv directory
        makePath((outdir + "/" + tsvdir).c_str());
    }

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
    tsv_reader ftr_tr((in_sgedir + "/" + sge_ftrf).c_str());
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
    notice("Processing SGE directory %s ...", in_sgedir.c_str());
    sge2_stream_reader ssr; // ((sgedir + "/" + bcdf).c_str(), (sgedir + "/" + ftrf).c_str(), (sgedir + "/" + mtxf).c_str());

    // load the position file if provided
    if (!posf.empty())
    {
        notice("Loading position file %s ...", posf.c_str());
        ssr.load_position_file(posf.c_str(), pos_colname_barcode.c_str(), pos_colname_x.c_str(), pos_colname_y.c_str(), pos_delim[0]);
    }
    notice("Opening the SGE directory %s ...", in_sgedir.c_str());

    // specify the input columns
    ssr.icol_bcd_strid = in_icol_bcd_barcode - 1;
    ssr.icol_bcd_px = in_icol_bcd_px - 1;
    ssr.icol_bcd_py = in_icol_bcd_py - 1;
    ssr.icol_ftr_id = in_icol_ftr_id - 1;
    ssr.icol_ftr_name = in_icol_ftr_name - 1;
    ssr.open((in_sgedir + "/" + sge_bcdf).c_str(), (in_sgedir + "/" + sge_ftrf).c_str(), (in_sgedir + "/" + sge_mtxf).c_str());

    htsFile* wh_tsv = out_tsv ? hts_open((outdir + "/" + tsvdir + "/" + tsv_mtxf).c_str(), "wz") : NULL;
    if ( ( out_tsv ) && ( wh_tsv == NULL ) ) 
        error("Failed to open %s/%s/%s for writing", outdir.c_str(), tsvdir.c_str(), tsv_mtxf.c_str());
    
    sge2_stream_writer* pssw = NULL;
    if ( out_sge ) {
        pssw = new sge2_stream_writer((outdir + "/" + sgedir + "/" + sge_bcdf).c_str(), 
            (outdir + "/" + sgedir + "/" + sge_ftrf).c_str(), 
            (outdir + "/" + sgedir + "/" + sge_mtxf).c_str());
    }

    // print the header
    if ( wh_tsv != NULL ) {
        if (print_feature_id)
            hprintf(wh_tsv, "%s\t%s\t%s\t%s\t%s\t%s\n", colname_x.c_str(), colname_y.c_str(),
                colname_gene_name.c_str(), colname_gene_id.c_str(), colname_count.c_str());
        else
            hprintf(wh_tsv, "%s\t%s\t%s\t%s\n", colname_x.c_str(), colname_y.c_str(),
                colname_gene_name.c_str(), colname_count.c_str());
    }

    // read each line of the matrix file
    double um_x = 0.0, um_y = 0.0;
    double xmax = -1e99, ymax = -1e99, xmin = 1e99, ymin = 1e99;
    uint64_t intervals = 1000000;
    uint64_t nlines = 0;
    while (ssr.read_mtx())
    {
        if (nlines % intervals == 0)
            notice("Processing %llu lines... from %s/%s", nlines, in_sgedir.c_str(), sge_mtxf.c_str());
        // remove lane and tile info and create global coordinates
        if (ssr.is_bcd_new)
        { // new barcode
            // calculate the new pixel coordinatex
            //notice("foo");
            um_x = ssr.cur_sbcd.px / units_per_um;
            um_y = ssr.cur_sbcd.py / units_per_um;

            // update the bounding box
            xmax = xmax > um_x ? xmax : um_x;
            ymax = ymax > um_y ? ymax : um_y;
            xmin = xmin < um_x ? xmin : um_x;
            ymin = ymin < um_y ? ymin : um_y;

            if ( out_sge ) {
                pssw->add_sbcd(ssr.cur_sbcd.strid.c_str(), um_x, um_y);
            }
        }
        // print individual transcripts
        if (ftr_passes[ssr.cur_iftr-1]) { // print only if the feature passes the filter
            if ( wh_tsv != NULL ) {
                if ( floor(um_x) != um_x || floor(um_y) != um_y ) { // non-integer coordinates, use the precision
                    hprintf(wh_tsv, "%.*f\t%.*f\t%s", precision_um, um_x, precision_um, um_y, ftr_names[ssr.cur_iftr-1].c_str());
                }
                else { // integer coordinates, ignore the precision
                    hprintf(wh_tsv, "%d\t%d\t%s", (int)floor(um_x), (int)floor(um_y), ftr_names[ssr.cur_iftr-1].c_str());
                }
                if (print_feature_id)
                    hprintf(wh_tsv, "\t%s", ftr_ids[ssr.cur_iftr-1].c_str());

                // print the count information
                hprintf(wh_tsv, "\t%llu\n", ssr.cur_cnts[in_icol_mtx-1]);
            }
        }
        ftr_cnts[ssr.cur_iftr-1] += ssr.cur_cnts[in_icol_mtx-1];
        if ( out_sge ) {
            pssw->add_mtx(ssr.cur_iftr, ssr.cur_cnts);
        }
        ++nlines;
    }
    if ( wh_tsv != NULL ) {
        hts_close(wh_tsv);
    }

    if ( out_tsv ) {
        htsFile* wh_ftr = hts_open((outdir + "/" + tsvdir + "/" + tsv_ftrf).c_str(), "wz");
        if ( wh_ftr == NULL )
            error("Failed to open %s/%s/%s for writing", outdir.c_str(), tsvdir.c_str(), tsv_ftrf.c_str());
        hprintf(wh_ftr, "%s\t%s\t%s\n", colname_gene_name.c_str(), colname_gene_id.c_str(), colname_count.c_str());
        for (int32_t i = 0; i < ftr_ids.size(); ++i)
        {
            if (ftr_passes[i] && ftr_cnts[i] > 0)
                hprintf(wh_ftr, "%s\t%s\t%llu\n", ftr_names[i].c_str(), ftr_ids[i].c_str(), ftr_cnts[i]);
        }
        hts_close(wh_ftr);

        htsFile* wh_minmax = hts_open((outdir + "/" + tsvdir + "/" + tsv_minmaxf).c_str(), "w");
        if ( wh_minmax == NULL )
            error("Failed to open %s/%s/%s for writing", outdir.c_str(), tsvdir.c_str(), tsv_minmaxf.c_str());
        hprintf(wh_minmax, "xmin\t%.*f\n", precision_um, xmin);
        hprintf(wh_minmax, "xmax\t%.*f\n", precision_um, xmax);
        hprintf(wh_minmax, "ymin\t%.*f\n", precision_um, ymin);
        hprintf(wh_minmax, "ymax\t%.*f\n", precision_um, ymax);
        hts_close(wh_minmax);
    }
    if ( out_sge ) {
        int32_t nftrs = ssr.load_features();
        pssw->ftr_cnts.resize(nftrs);
        for (int32_t j = 0; j < nftrs; ++j)
        {
            pssw->write_ftr(ssr.ftrs[j]->id.c_str(), ssr.ftrs[j]->name.c_str(), (uint64_t)(j + 1), pssw->ftr_cnts[j]);
        }
        pssw->close();
        delete pssw;
    }

    ssr.close();

    notice("Analysis finished");

    return 0;
}
