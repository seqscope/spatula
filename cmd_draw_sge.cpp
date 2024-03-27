#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge.h"
#include <cmath>
#include <ctime>
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
        if (s_color[0] != '#' || strlen(s_color) != 7)
            return false;
        sscanf(s_color + 1, "%02x%02x%02x", &rgb[0], &rgb[1], &rgb[2]);
        return true;
    }
};
typedef struct _color_gene_unit_t color_gene_unit_t;

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
    double coord_per_pixel = 1000.0; // 1 pixel = 1 coordinate unit
    std::string outf;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input files", NULL)
    LONG_STRING_PARAM("manifest", &manifestf, "Manifest file from combine-sbcds. xmin/xmax/ymin/ymax will be automatically detected")
    LONG_STRING_PARAM("sge", &sgedir, "SGE directory")
    LONG_STRING_PARAM("bcd", &bcdf, "Barcode file name")
    LONG_STRING_PARAM("ftr", &ftrf, "Feature file name")
    LONG_STRING_PARAM("mtx", &mtxf, "Matrix file name")

    LONG_PARAM_GROUP("Genes to visualize", NULL)
    LONG_MULTI_STRING_PARAM("color-gene", &color_genes, "[color_code]:[gene1],[gene2],... as a visualization unit. Adding :[idx] at the end is optional")
    LONG_MULTI_STRING_PARAM("color-list", &color_lists, "[color_code]:[list_file] as a visualization unit")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_DOUBLE_PARAM("coord-per-pixel", &coord_per_pixel, "Number of coordinate units per pixel")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file name")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (manifestf.empty() || outf.empty() || sgedir.empty())
        error("--manifest, --sge and --out must be specified");

    notice("Analysis started");

    // read the manifest file and determine the xmin/xmax/ymin/ymax
    uint64_t xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    dataframe_t df(manifestf.c_str());
    if (df.nrows == 0)
        error("Empty dataframe %s", manifestf.c_str());
    for (int32_t i = 0; df.nrows; ++i)
    {
        // parse xmin/xmax/ymin/ymax
        int32_t i_xmin = df.get_colidx("xmin");
        int32_t i_xmax = df.get_colidx("xmax");
        int32_t i_ymin = df.get_colidx("ymin");
        int32_t i_ymax = df.get_colidx("ymax");
        if (i_xmin < 0 || i_xmax < 0 || i_ymin < 0 || i_ymax < 0)
            error("Cannot find the columns 'xmin', 'xmax', 'ymin', 'ymax' in %s", manifestf.c_str());
        xmin = df.get_uint64_elem(i, i_xmin);
        xmax = df.get_uint64_elem(i, i_xmax);
        ymin = df.get_uint64_elem(i, i_ymin);
        ymax = df.get_uint64_elem(i, i_ymax);
        break;
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
            for (int32_t j = 0; j < v2.size(); ++j)
            {
                auto range = ftr_name2idx.equal_range(v2[j]);
                std::vector<uint32_t> igenes;
                if (range.first == range.second)
                { // gene name is not recognized
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
                for (int32_t j = 0; j < igenes.size(); ++j)
                {
                    cgu.igene2col[igenes[j]] = idx;
                    igene2icgus[igenes[j]].push_back(i);
                }
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
            std::vector<uint32_t> igenes;
            if (icol_id >= 0) // gene id is available
            {
                const std::string& gid = df.get_str_elem(j, icol_id);
                if (ftr_id2idx.find(gid) == ftr_id2idx.end())
                {
                    error("Cannot find gene ID %s from the SGE matrix", gid.c_str());
                }
                igenes.push_back(ftr_id2idx[gid]);
            }
            else
            {
                const std::string& gname = df.get_str_elem(j, icol_name);
                auto range = ftr_name2idx.equal_range(gname);
                if (range.first == range.second)
                { // gene name is not recognized
                    error("Cannot recognize gene name %s from the SGE matrix", gname.c_str());
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
                igene2icgus[igenes[j]].push_back(i + color_genes.size());
            }
        }
        color_gene_units.push_back(cgu);
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
        if (igene2icgus.find(ssr.cur_iftr) != igene2icgus.end())
        {
            for (int32_t i = 0; i < igene2icgus[ssr.cur_iftr].size(); ++i)
            {
                int32_t icgu = igene2icgus[ssr.cur_iftr][i]; // color gene unit index
                int32_t cnt = (int32_t)ssr.cur_cnts[color_gene_units[icgu].igene2col[ssr.cur_iftr]];
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

    notice("Writing the image to %s", outf.c_str());
    image.save_png(outf.c_str());

    return 0;
}
