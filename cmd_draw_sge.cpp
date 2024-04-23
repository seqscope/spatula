#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#define cimg_display 0 // remove the need for X11 library
#include "cimg/CImg.h"

struct _color_gene_unit_t
{
    int32_t rgb[3];
    std::map<uint32_t, int32_t> igene2col; // gene index to column index
    _color_gene_unit_t()
    {
        rgb[0] = rgb[1] = rgb[2] = 0;
    }
    bool set_rgb(const char *s_color)
    {
        if (s_color[0] == '#' && strlen(s_color) == 7) {
            sscanf(s_color + 1, "%02x%02x%02x", &rgb[0], &rgb[1], &rgb[2]);
            return true;
        }
        else if ( strlen(s_color) == 6 ) {
            sscanf(s_color, "%02x%02x%02x", &rgb[0], &rgb[1], &rgb[2]);
            return true;
        }
        else {
            error("Invalid color code %s", s_color);
            return false;
        }
    }
};
typedef struct _color_gene_unit_t color_gene_unit_t;

uint32_t get_quantile(std::map<uint64_t, uint32_t> &hist, double quantile)
{
    uint64_t total = 0;
    std::map<uint32_t, uint64_t> cnts;
    for (auto it = hist.begin(); it != hist.end(); ++it)
    {
        ++total;// += it->second;
        ++cnts[it->second];
    }
    uint64_t sum = 0;
    for (auto it = cnts.begin(); it != cnts.end(); ++it)
    {
        sum += it->second;
        if (sum >= total * quantile)
            return it->first;
    }
    warning("Cannot find the quantile %f", quantile);
    return 1;
}

void adjust_rgb_by_max(int32_t* rgb, uint32_t maxval, int32_t cap = 255) {
    int32_t maxr = rgb[0] * maxval;
    int32_t maxg = rgb[1] * maxval;
    int32_t maxb = rgb[2] * maxval;
    int32_t maxc = std::max(std::max(maxr, maxg), maxb);
    double ratio = (double)cap / (maxc > 0 ? maxc : 1);
    rgb[0] = (int32_t)(rgb[0] * ratio);
    rgb[1] = (int32_t)(rgb[1] * ratio);
    rgb[2] = (int32_t)(rgb[2] * ratio);
}

/////////////////////////////////////////////////////////////////////////
// draw-xy : Draw the single-color image of points in 2D space
////////////////////////////////////////////////////////////////////////
int32_t cmdDrawSGE(int32_t argc, char **argv)
{
    std::string manifestf;
    std::string sgedir;
    std::string bcdf("barcodes.tsv.gz");
    std::string ftrf("features.tsv.gz");
    std::string mtxf("matrix.mtx.gz");
    std::vector<std::string> color_genes;
    std::vector<std::string> color_lists;
    std::vector<std::string> color_regexes;
    double coord_per_pixel = 1000.0; // 1 pixel = 1 coordinate unit
    bool auto_adjust_intensity = false;
    double auto_adjust_quantile = 0.99;
    std::string outf;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input files", NULL)
    LONG_STRING_PARAM("minmax", &manifestf, "Bounding box information. Expects xmin/xmax/ymin/ymax (tall or wide format)")
    LONG_STRING_PARAM("sge", &sgedir, "SGE directory")
    LONG_STRING_PARAM("bcd", &bcdf, "Barcode file name")
    LONG_STRING_PARAM("ftr", &ftrf, "Feature file name")
    LONG_STRING_PARAM("mtx", &mtxf, "Matrix file name")

    LONG_PARAM_GROUP("Genes to visualize", NULL)
    LONG_MULTI_STRING_PARAM("color-gene", &color_genes, "[color_code]:[gene1],[gene2],... as a visualization unit. Adding :[idx] at the end is optional")
    LONG_MULTI_STRING_PARAM("color-regex", &color_regexes, "[color_code]:[regex](:[idx]) as a visualization unit. [regex] is a regulat expression. Adding :[idx] at the end is optional")
    LONG_MULTI_STRING_PARAM("color-list", &color_lists, "[color_code]:[list_file] as a visualization unit")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_DOUBLE_PARAM("coord-per-pixel", &coord_per_pixel, "Number of coordinate units per pixel")
    LONG_PARAM("auto-adjust", &auto_adjust_intensity, "Automatically adjust the intensity of the color based on the maximum count")
    LONG_DOUBLE_PARAM("adjust-quantile", &auto_adjust_quantile, "Quantile of pixel to use for auto-adjustment among non-zero pixels")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file name")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (outf.empty() || sgedir.empty())
        error("--sge and --out must be specified");

    notice("Analysis started");

    // read the manifest file and determine the xmin/xmax/ymin/ymax
    uint64_t xmin = UINT64_MAX, xmax = 0, ymin = 0, ymax = UINT64_MAX;
    if ( !manifestf.empty() ) {
        read_minmax(manifestf.c_str(), xmin, xmax, ymin, ymax);
    }
    else if ( !auto_adjust_intensity ) {
        error("Cannot determine the bounding box. Please specify the manifest file or turn on the --auto-adjust option");
    }

    // parse the gene list information
    tsv_reader ftr_tr((sgedir + "/" + ftrf).c_str());
    std::vector<std::string> ftr_ids;
    std::vector<std::string> ftr_names;
    std::multimap<std::string, uint32_t> ftr_name2idx;
    std::map<std::string, uint32_t> ftr_id2idx;
    while (ftr_tr.read_line() > 0)
    {
        ftr_ids.push_back(ftr_tr.str_field_at(0));
        ftr_names.push_back(ftr_tr.str_field_at(1));
    }
    for (int32_t i = 0; i < ftr_ids.size(); ++i)
    {
        ftr_name2idx.insert(std::pair<std::string, uint32_t>(ftr_names[i], i));
        ftr_id2idx[ftr_ids[i]] = i;
    }

    std::vector<color_gene_unit_t> color_gene_units;
    std::map<int32_t, std::vector<int32_t>> igene2icgus; // gene index to color gene unit index

    // parse color-genes arguments first
    for (int32_t i = 0; i < color_genes.size(); ++i)
    {
        color_gene_unit_t cgu;
        std::vector<std::string> v1, v2;
        split(v1, ":", color_genes[i]);
        cgu.set_rgb(v1[0].c_str());
        if (v1.size() < 2)
            error("Invalid color-gene argument: %s", color_genes[i].c_str());
        int32_t idx = v1.size() > 2 ? atoi(v1[2].c_str()) - 1 : 0;
        if (v1[1].compare("_all_") == 0)
        { // if "_all_" is specified include all genes
            for (int32_t j = 0; j < (int32_t)ftr_ids.size(); ++j)
            {
                cgu.igene2col[j] = idx;
                igene2icgus[j].push_back(i);
            }
        }
        else
        { // otherwise, parse the gene names
            split(v2, ",", v1[1]);
            int32_t ngenes = 0;
            for (int32_t j = 0; j < v2.size(); ++j)
            {
                auto range = ftr_name2idx.equal_range(v2[j]);
                std::vector<uint32_t> igenes;
                if (range.first == range.second)
                {   // gene name is not recognized, spit out an error
                    if (ftr_id2idx.find(v2[j]) == ftr_id2idx.end())
                    {
                        error("Cannot recognize %s as gene name or ID", v2[j].c_str());
                    }
                    igenes.push_back(ftr_id2idx[v2[j]]);
            
                }
                else
                {
                    for (auto it = range.first; it != range.second; ++it)
                    {
                        igenes.push_back(it->second);
                    }
                    if (igenes.size() > 1)
                        warning("Gene name %s have %d matching gene ID", v2[j].c_str(), (int32_t)igenes.size());
                }
                ngenes += igenes.size();
                for (int32_t j = 0; j < igenes.size(); ++j)
                {
                    cgu.igene2col[igenes[j]] = idx;
                    igene2icgus[igenes[j]].push_back(i);
                }
            }
            if (ngenes == 0)
                error("No gene is found from the gene list: %s", v1[1].c_str());
            notice("Identified %d genes from the gene list: %s", ngenes, v1[1].c_str());
        }
        color_gene_units.push_back(cgu);
    }

    // parse color-regex arguments second
    for (int32_t i = 0; i < color_regexes.size(); ++i)
    {
        color_gene_unit_t cgu;
        std::vector<std::string> v1, v2;
        split(v1, ":", color_regexes[i]);
        cgu.set_rgb(v1[0].c_str());
        if (v1.size() < 2)
            error("Invalid color-regex argument: %s", color_regexes[i].c_str());
        int32_t idx = v1.size() > 2 ? atoi(v1[2].c_str()) - 1 : 0;

        // scan each gene name and check if it matches the regular expression
        std::regex re(v1[1]);
        for (int32_t j = 0; j < (int32_t)ftr_ids.size(); ++j)
        {
            if ( std::regex_match(ftr_names[j], re) ) {
                cgu.igene2col[j] = idx;
                igene2icgus[j].push_back(i + color_genes.size());
            }
        }
        color_gene_units.push_back(cgu);
    }

    // parse color-genes list
    // possible headers are [id] [name] [idx]
    for (int32_t i = 0; i < color_lists.size(); ++i)
    {
        color_gene_unit_t cgu;
        std::vector<std::string> v1;
        split(v1, ":", color_lists[i]);
        if (v1.size() < 2)
            error("Invalid color-list argument: %s", color_lists[i].c_str());
        cgu.set_rgb(v1[0].c_str());
        dataframe_t df(v1[1].c_str());
        int32_t default_idx = v1.size() > 2 ? atoi(v1[2].c_str()) - 1 : 0;
        int32_t icol_id = df.get_colidx("id");
        int32_t icol_name = df.get_colidx("name");
        int32_t icol_idx = df.get_colidx("idx");
        if ((icol_id < 0) && (icol_name < 0))
        {
            error("Cannot find the column 'id' or 'name' in %s", v1[1].c_str());
        }
        for (int32_t j = 0; j < df.nrows; ++j)
        {
            bool warning_no_match = false;
            std::vector<uint32_t> igenes;
            if (icol_id >= 0) // gene id is available
            {
                const std::string& gid = df.get_str_elem(j, icol_id);
                // if gene name is NOT recognized, then pass. Spit error only when no gene is recognized
                if (ftr_id2idx.find(gid) == ftr_id2idx.end())
                {
                    //do nothing, just skip this gene
                    //error("Cannot find gene ID %s from the SGE matrix", gid.c_str());
                    warning_no_match = true;
                }
                else {
                    igenes.push_back(ftr_id2idx[gid]);
                }
            }
            else
            {
                const std::string& gname = df.get_str_elem(j, icol_name);
                auto range = ftr_name2idx.equal_range(gname);
                if (range.first == range.second)
                {   // gene name is not recognized
                    warning_no_match = true;
                    // error("Cannot recognize gene name %s from the SGE matrix", gname.c_str());
                }
                for (auto it = range.first; it != range.second; ++it)
                {
                    igenes.push_back(it->second);
                }
                if (igenes.size() > 1)
                    warning("Gene name %s have %d matching gene ID", gname.c_str(), (int32_t)igenes.size());
            }
            int32_t idx = icol_idx >= 0 ? df.get_int_elem(j, icol_idx) - 1 : default_idx;
            for (int32_t j = 0; j < igenes.size(); ++j)
            {
                cgu.igene2col[igenes[j]] = idx;
                igene2icgus[igenes[j]].push_back(i + color_genes.size() + color_regexes.size());
            }
        }
        color_gene_units.push_back(cgu);
    }

    if (color_gene_units.empty())
        error("No color-gene units are specified. Please use at least one --color-gene, --color-regex or --color-list option");

    notice("Parsed color-gene units:");
    for (int32_t i = 0; i < color_gene_units.size(); ++i)
    {
        notice("Color gene unit %d: RGB(%x,%x,%x)", i, color_gene_units[i].rgb[0], color_gene_units[i].rgb[1], color_gene_units[i].rgb[2]);
        if ( color_gene_units[i].igene2col.size() < 10 ) {
            for (auto it = color_gene_units[i].igene2col.begin(); it != color_gene_units[i].igene2col.end(); ++it)
            {
                notice("  Gene %d : %s (%s) -> Color index %d", it->first, ftr_ids[it->first].c_str(), ftr_names[it->first].c_str(), it->second);
            }
        }
        else {
            notice("  ... %d genes (too many to show individual genes)", (int32_t)color_gene_units[i].igene2col.size());
        }
    }

    // if auto-adjust is ON, then create histogram for individual colors
    if (auto_adjust_intensity)
    {
        notice("Auto-adjusting the intensity of the colors... Pre-scanning the SGE matrix...");
        std::vector< std::map<uint64_t, uint32_t> > color_hist(color_gene_units.size());
        sge_stream_reader ssr((sgedir + "/" + bcdf).c_str(), (sgedir + "/" + ftrf).c_str(), (sgedir + "/" + mtxf).c_str());

        // read each line of the matrix file
        int32_t px = -1, py = -1;
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
                px = ssr.cur_sbcd.px / coord_per_pixel; // do not subtract xmin/ymin
                py = ssr.cur_sbcd.py / coord_per_pixel;
                // adjust bounding box
                if ( auto_adjust_intensity && manifestf.empty() ) {
                    if (ssr.cur_sbcd.px < xmin) xmin = ssr.cur_sbcd.px;
                    if (ssr.cur_sbcd.px > xmax) xmax = ssr.cur_sbcd.px;
                    if (ssr.cur_sbcd.py < ymin) ymin = ssr.cur_sbcd.py;
                    if (ssr.cur_sbcd.py > ymax) ymax = ssr.cur_sbcd.py;
                }
            }
            // check if the gene is relevant
            if (igene2icgus.find(ssr.cur_iftr-1) != igene2icgus.end())
            {
                for (int32_t i = 0; i < igene2icgus[ssr.cur_iftr-1].size(); ++i)
                {
                    int32_t icgu = igene2icgus[ssr.cur_iftr-1][i]; // color gene unit index
                    int32_t cnt = (int32_t)ssr.cur_cnts[color_gene_units[icgu].igene2col[ssr.cur_iftr-1]];
                    if (cnt > 0) {
                        uint64_t pxy = ((uint64_t)px << 32) | py;
                        color_hist[icgu][pxy] += cnt;
                    }
                }
            }
            ++nlines;
        }
        ssr.close();

        // based on the histogram, update the color profile,
        for (int32_t i = 0; i < color_gene_units.size(); ++i)
        {
            int32_t rgb[3] = {color_gene_units[i].rgb[0], color_gene_units[i].rgb[1], color_gene_units[i].rgb[2]};
            uint32_t maxval = get_quantile(color_hist[i], auto_adjust_quantile);
            adjust_rgb_by_max(rgb, maxval, 255);
            notice("Color gene unit %d: RGB(%x,%x,%x) -> Adjusted RGB(%x,%x,%x) with maxval = %d", i, color_gene_units[i].rgb[0], color_gene_units[i].rgb[1], color_gene_units[i].rgb[2], rgb[0], rgb[1], rgb[2], maxval);
            color_gene_units[i].rgb[0] = rgb[0];
            color_gene_units[i].rgb[1] = rgb[1];
            color_gene_units[i].rgb[2] = rgb[2];
        }
    }

    int32_t height = (int32_t)(ceil((double)(ymax - ymin + 1) / coord_per_pixel));
    int32_t width = (int32_t)(ceil((double)(xmax - xmin + 1) / coord_per_pixel));
    cimg_library::CImg<unsigned char> image(width, height, 1, 3, 0);

    // read the SGE matrix
    notice("Processing SGE directory %s ...", sgedir.c_str());
    sge_stream_reader ssr((sgedir + "/" + bcdf).c_str(), (sgedir + "/" + ftrf).c_str(), (sgedir + "/" + mtxf).c_str());
    // read each line of the matrix file
    int32_t px = -1, py = -1;
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
            px = (ssr.cur_sbcd.px - xmin) / coord_per_pixel;
            py = (ssr.cur_sbcd.py - ymin) / coord_per_pixel;
        }
        // check if the gene is relevant
        if (igene2icgus.find(ssr.cur_iftr-1) != igene2icgus.end())
        {
            //notice("%llu %llu", ssr.cur_sbcd.nid, ssr.cur_iftr);
            for (int32_t i = 0; i < igene2icgus[ssr.cur_iftr-1].size(); ++i)
            {
                int32_t icgu = igene2icgus[ssr.cur_iftr-1][i]; // color gene unit index
                int32_t cnt = (int32_t)ssr.cur_cnts[color_gene_units[icgu].igene2col[ssr.cur_iftr-1]];
                if (cnt > 0)
                {
                    // update the color at the pixel
                    for (int32_t j = 0; j < 3; ++j)
                    { //
                        int32_t c = image(px, py, j);
                        if (c + color_gene_units[icgu].rgb[j] > 255)
                            c = 255;
                        else
                            c += color_gene_units[icgu].rgb[j];
                        image(px, py, j) = c;
                    }
                }
            }
        }
        ++nlines;
    }
    ssr.close();

    notice("Writing the image to %s", outf.c_str());
    image.save_png(outf.c_str());

    return 0;
}
