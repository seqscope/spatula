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
    //std::string outdir;
    std::string out_sgedir;
    std::string out_tsvdir;
    //bool out_sge = false;
    //bool out_tsv = false;
    std::string posf;   // external position files

    // SGE and TSV file names
    std::string sgedir("sge");
    std::string tsvdir("tsv");
    std::string sge_bcdf("barcodes.tsv.gz");
    std::string sge_ftrf("features.tsv.gz");
    std::string sge_mtxf("matrix.mtx.gz");
    std::string tsv_mtxf("transcripts.unsorted.tsv.gz");
    std::string tsv_ftrf("features.tsv.gz");
    std::string tsv_minmaxf("minmax.tsv");

    // input column indices
    int32_t in_icol_bcd_barcode = 1;
    int32_t in_icol_bcd_px = 6; // 1-based column index in the barcode file to use as the X coordinate
    int32_t in_icol_bcd_py = 7; // 1-based column index in the barcode file to use as the Y coordinate
    int32_t in_icol_ftr_id = 1;
    int32_t in_icol_ftr_name = 2;
    int32_t in_icol_mtx_thres = -1; 
    int32_t mtx_thres = 0; // threshold for the matrix file to be considered as a valid count
    // int32_t in_icol_mtx = 1;
    std::string str_icols_mtx("1,2,3,4,5");
    std::vector<int32_t> v_icols_mtx;

    // output column names for tsv files
    std::string colname_gene_name("gene");
    std::string colname_gene_id("gene_id");
    //std::string colname_count("Count");
    std::string colnames_count("gn,gt,spl,unspl,ambig");
    std::vector<std::string> v_colnames_count;
    std::string colname_x("X");
    std::string colname_y("Y");
    bool print_feature_id = false;
    bool allow_duplicate_gene_names = false;

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
    LONG_STRING_PARAM("out-sge", &out_sgedir, "Path of SGE output directory")
    LONG_STRING_PARAM("out-tsv", &out_tsvdir, "Path of TSV output directory")

    LONG_PARAM_GROUP("File names for SGE/TSV input/output", NULL)
    LONG_STRING_PARAM("sge-bcd", &sge_bcdf, "Barcode file name in SGE directory")
    LONG_STRING_PARAM("sge-ftr", &sge_ftrf, "Feature file name in SGE directory")
    LONG_STRING_PARAM("sge-mtx", &sge_mtxf, "Matrix file name in SGE directory")
    LONG_STRING_PARAM("tsv-mtx", &tsv_mtxf, "Transcript file name in TSV directory")
    LONG_STRING_PARAM("tsv-ftr", &tsv_ftrf, "Feature file name in TSV directory")
    LONG_STRING_PARAM("tsv-minmax", &tsv_minmaxf, "Minmax file name in TSV directory")

    LONG_PARAM_GROUP("Expected column index in SGE input", NULL)
    LONG_STRING_PARAM("icols-mtx", &str_icols_mtx, "Comma-separated 1-based column indices use as the count")
    LONG_INT_PARAM("icol-bcd-barcode", &in_icol_bcd_barcode, "1-based column index of barcode in the barcode file")
    LONG_INT_PARAM("icol-bcd-x", &in_icol_bcd_px, "1-based column index of x coordinate in the barcode file")
    LONG_INT_PARAM("icol-bcd-y", &in_icol_bcd_py, "1-based column index of y coordinate in the barcode file")
    LONG_INT_PARAM("icol-ftr-id", &in_icol_ftr_id, "1-based column index of feature ID in the barcode file")
    LONG_INT_PARAM("icol-ftr-name", &in_icol_ftr_name, "1-based column index of feature name in the barcode file")

    LONG_PARAM_GROUP("Thrsholding options", NULL)
    LONG_INT_PARAM("icol-thres", &in_icol_mtx_thres, "1-based column index of the threshold in the matrix file (default: -1)")    
    LONG_INT_PARAM("mtx-thres", &mtx_thres, "Threshold for the matrix file to be considered as a valid count (default: 0)")

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
    LONG_PARAM("allow-duplicate-gene-names", &allow_duplicate_gene_names, "Allow duplicate gene names in the output file")
    LONG_STRING_PARAM("colname-feature-name", &colname_gene_name, "Column name for feature/gene name")
    LONG_STRING_PARAM("colname-feature-id", &colname_gene_id, "Column name for feature/gene ID")
    LONG_STRING_PARAM("colnames-count", &colnames_count, "Comma-separate column names for Count")
    LONG_STRING_PARAM("colname-x", &colname_x, "Column name for X")
    LONG_STRING_PARAM("colname-y", &colname_y, "Column name for Y")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( in_sgedir.empty() || ( out_tsvdir.empty() && out_sgedir.empty() ) )
        error("Input/output directory --in-sge and --out-dir must be specified");

    if ( out_sgedir.empty() && out_tsvdir.empty() )
        error("At least one of --out-sge or --out-tsv must be specified");

    // strip the trailing slash
    if ( ( !in_sgedir.empty() ) && ( in_sgedir.back() == '/' ) )
        in_sgedir.pop_back();

    if ( ( !out_sgedir.empty()) && ( out_sgedir.back() == '/' ) )
        out_sgedir.pop_back();

    if ( ( !out_tsvdir.empty() ) && ( out_tsvdir.back() == '/' ) )
        out_tsvdir.pop_back();

    if ( !out_sgedir.empty() ) { // create sge directory
        if ( isDirExist(out_sgedir) ) {
            notice("Output directory %s already exists.", out_sgedir.c_str());
        }
        else {   
            notice("Creating the directory %s....", out_sgedir.c_str());
            makePath(out_sgedir);
        }
    }
    if ( !out_tsvdir.empty() ) {   // create tsv directory
        if ( isDirExist(out_tsvdir) ) {
            notice("Output directory %s already exists.", out_tsvdir.c_str());
        }
        else {   
            notice("Creating the directory %s....", out_tsvdir.c_str());
            makePath(out_tsvdir);
        }
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

    // parse the column indices
    std::vector<std::string> v_str_icols_mtx;
    split(v_str_icols_mtx, ",", str_icols_mtx);
    split(v_colnames_count, ",", colnames_count);
    if ( v_str_icols_mtx.size() != v_colnames_count.size() )
        error("--icols-mtx [%s] and --colnames-count [%s] has different number of entries", str_icols_mtx.c_str(), colnames_count.c_str());
    for (int32_t i = 0; i < v_str_icols_mtx.size(); ++i)
    {
        v_icols_mtx.push_back(atoi(v_str_icols_mtx[i].c_str())-1);
    }
    int32_t n_colnames = (int32_t)v_colnames_count.size();

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
    std::vector< std::vector<uint64_t> > ftr_cnts;
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
        ftr_cnts.resize(ftr_cnts.size() + 1);
        ftr_cnts.back().resize(n_colnames, 0);
    }

    // ensure that feature IDs are unique
    std::set<std::string> ftr_id_set;
    for (int32_t i = 0; i < ftr_ids.size(); ++i)
    {
        if (ftr_id_set.find(ftr_ids[i]) != ftr_id_set.end())
            error("Feature ID %s is not unique", ftr_ids[i].c_str());
        ftr_id_set.insert(ftr_ids[i]);
    }

    if ( !allow_duplicate_gene_names ) {
        // force the feature name to be unique
        std::map<std::string, std::vector<int32_t> > ftr_name_indices;
        for (int32_t i = 0; i < ftr_names.size(); ++i)
        {
            ftr_name_indices[ftr_names[i]].push_back(i);
        }

        // resolve duplicate feature names
        for (std::map<std::string, std::vector<int32_t> >::iterator it = ftr_name_indices.begin(); it != ftr_name_indices.end(); ++it)
        {
            if (it->second.size() > 1)
            {
                for (int32_t i = 0; i < it->second.size(); ++i)
                {
                    notice("Feature name %s (ID : %s) is not unique. Resolving the conflict by appending the suffix v%d", it->first.c_str(), ftr_ids[it->second[i]].c_str(), i + 1);
                    ftr_names[it->second[i]] = ftr_names[it->second[i]] + "_v" + std::to_string(i + 1);
                    //ftr_names[it->second[i]] = ftr_names[it->second[i]] + "_" + ftr_ids[it->second[i]];
                }
            }
        }
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

    htsFile* wh_tsv = out_tsvdir.empty() ? NULL : hts_open((out_tsvdir + "/" + tsv_mtxf).c_str(), "wz");
    if ( ( !out_tsvdir.empty() ) && ( wh_tsv == NULL ) ) 
        error("Failed to open %s/%s/ for writing", out_tsvdir.c_str(), tsv_mtxf.c_str());
    
    sge2_stream_writer* pssw = NULL;
    if ( !out_sgedir.empty() ) {
        pssw = new sge2_stream_writer((out_sgedir + "/" + sge_bcdf).c_str(), 
            (out_sgedir + "/" + sge_ftrf).c_str(), 
            (out_sgedir + "/" + sge_mtxf).c_str());
    }

    // print the header
    if ( wh_tsv != NULL ) {
        hprintf(wh_tsv, "%s\t%s\t%s", colname_x.c_str(), colname_y.c_str(), colname_gene_name.c_str());
        if (print_feature_id)
            hprintf(wh_tsv, "\t%s", colname_gene_id.c_str());
        for (int32_t i = 0; i < n_colnames; ++i)
        {
            hprintf(wh_tsv, "\t%s", v_colnames_count[i].c_str());
        }
        hprintf(wh_tsv, "\n");
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

            if ( pssw != NULL ) {
                pssw->add_sbcd(ssr.cur_sbcd.strid.c_str(), um_x, um_y);
            }
        }
        // print individual transcripts
        if (ftr_passes[ssr.cur_iftr-1]) { // print only if the feature passes the filter
            uint64_t sum_cur_cnts = 0;
            for(int32_t i=0; i < n_colnames; ++i) {
                sum_cur_cnts += ssr.cur_cnts[v_icols_mtx[i]];
            }
            if ( sum_cur_cnts == 0 ) {
                // skip the feature if the count is zero
                continue;
            }

            // if the threshold is specified, skip the feature if the count is below the specified threshold
            if ( in_icol_mtx_thres > 0 ) {
                if ( ssr.cur_cnts[in_icol_mtx_thres-1] < mtx_thres ) {
                    continue;
                }
            }                   

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
                //hprintf(wh_tsv, "\t%llu\n", ssr.cur_cnts[in_icol_mtx-1]);
                if ( n_colnames > ssr.cur_cnts.size() )
                    error("Number of columns in the count %d file is less than the number of columns in the matrix file %zu", n_colnames, ssr.cur_cnts.size());

                for (int32_t i = 0; i < n_colnames; ++i)
                {
                    hprintf(wh_tsv, "\t%llu", ssr.cur_cnts[v_icols_mtx[i]]);
                }
                hprintf(wh_tsv, "\n");
            }

            //ftr_cnts[ssr.cur_iftr-1] += ssr.cur_cnts[in_icol_mtx-1];
            for (int32_t i = 0; i < n_colnames; ++i)
            {
                ftr_cnts[ssr.cur_iftr-1][i] += ssr.cur_cnts[v_icols_mtx[i]];
            }
            if ( pssw != NULL) {
                pssw->add_mtx(ssr.cur_iftr, ssr.cur_cnts, v_icols_mtx);
            }
        }
        ++nlines;
    }
    if ( wh_tsv != NULL ) {
        hts_close(wh_tsv);
    }

    if ( !out_tsvdir.empty() ) {
        htsFile* wh_ftr = hts_open((out_tsvdir + "/" + tsv_ftrf).c_str(), "wz");
        if ( wh_ftr == NULL )
            error("Failed to open %s/%s for writing", out_tsvdir.c_str(), tsv_ftrf.c_str());

        hprintf(wh_ftr, "%s\t%s", colname_gene_name.c_str(), colname_gene_id.c_str());
        for (int32_t i = 0; i < n_colnames; ++i)
        {
            hprintf(wh_ftr, "\t%s", v_colnames_count[i].c_str());
        }
        hprintf(wh_ftr, "\n");

        for (int32_t i = 0; i < ftr_ids.size(); ++i)
        {
            bool print_gene = false;
            if ( ftr_passes[i] ) {
                for(int32_t j = 0; j < n_colnames; ++j) {
                    if (ftr_cnts[i][j] > 0) {
                        print_gene = true;
                        break;
                    }
                }
            }
            if (print_gene) {
                hprintf(wh_ftr, "%s\t%s", ftr_names[i].c_str(), ftr_ids[i].c_str());
                for(int32_t j = 0; j < n_colnames; ++j) {
                    hprintf(wh_ftr, "\t%llu", ftr_cnts[i][j]);
                }
                hprintf(wh_ftr, "\n");
            }
        }
        hts_close(wh_ftr);

        htsFile* wh_minmax = hts_open((out_tsvdir + "/" + tsv_minmaxf).c_str(), "w");
        if ( wh_minmax == NULL )
            error("Failed to open %s/%s for writing", out_tsvdir.c_str(), tsv_minmaxf.c_str());
        hprintf(wh_minmax, "xmin\t%.*f\n", precision_um, xmin);
        hprintf(wh_minmax, "xmax\t%.*f\n", precision_um, xmax);
        hprintf(wh_minmax, "ymin\t%.*f\n", precision_um, ymin);
        hprintf(wh_minmax, "ymax\t%.*f\n", precision_um, ymax);
        hts_close(wh_minmax);
    }
    if ( pssw != NULL) {
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
