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

uint32_t hex2rgba(const char *hex) {
    if (hex[0] == '#') {
        hex++;
    }
    if (strlen(hex) == 6) {
        int32_t r, g, b;
        sscanf(hex, "%02x%02x%02x", &r, &g, &b);
        return (r << 24) | (g << 16) | (b << 8) | 0xFF; // alpha is set to 255
    } else if (strlen(hex) == 8) {
        int32_t r, g, b, a;
        sscanf(hex, "%02x%02x%02x%02x", &r, &g, &b, &a);
        return (r << 24) | (g << 16) | (b << 8) | a; // use provided alpha
    } else {
        error("Invalid hex color code: %s", hex);
        return 0; // should not reach here
    }
}

/////////////////////////////////////////////////////////////////////////
// draw-tsv : Draw the single-color image of points in 2D space
////////////////////////////////////////////////////////////////////////
int32_t cmdDrawTSV(int32_t argc, char **argv)
{
    std::string minmaxf; // file containing xmin/xmax/ymin/ymax
    std::string tsvf;
    std::string cmapf; // gene name and RGB color
    std::string tsv_colname_X("X");
    std::string tsv_colname_Y("Y");
    std::string tsv_colname_gene("gene");
    std::string tsv_colname_count("count");
    std::string cmap_colname_gene("gene");
    std::string cmap_colname_hex("hex");
    double coord_per_pixel = 1.0; // 1 pixel = 1 coordinate unit
    double intensity_per_count = 0.1; // 1 count = 0.1 intensity in RGB color space
    bool skip_tsv_header = false; // skip the header in the TSV file
    std::string outf;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input files", NULL)
    LONG_STRING_PARAM("minmax", &minmaxf, "Bounding box information. Expects xmin/xmax/ymin/ymax")
    LONG_STRING_PARAM("tsv", &tsvf, "Input TSV file")
    LONG_STRING_PARAM("cmap", &cmapf, "Color map file containing gene name and RGB color")

    LONG_PARAM_GROUP("Column names", NULL)
    LONG_STRING_PARAM("tsv-colname-x", &tsv_colname_X, "Column name for the X coordinate in the input TSV file")
    LONG_STRING_PARAM("tsv-colname-y", &tsv_colname_Y, "Column name for the Y coordinate in the input TSV file")
    LONG_STRING_PARAM("tsv-colname-gene", &tsv_colname_gene, "Column name for the gene in the input TSV file")
    LONG_STRING_PARAM("tsv-colname-count", &tsv_colname_count, "Column name for the count in the input TSV file")
    LONG_STRING_PARAM("cmap-colname-gene", &cmap_colname_gene, "Column name for the gene in the color map file")
    LONG_STRING_PARAM("cmap-colname-hex", &cmap_colname_hex, "Column name for the hex color in the color map file")
    LONG_PARAM("skip-tsv-header", &skip_tsv_header, "Skip the header in the input TSV file")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_DOUBLE_PARAM("coord-per-pixel", &coord_per_pixel, "Number of coordinate units per pixel")
    LONG_DOUBLE_PARAM("intensity-per-count", &intensity_per_count, "Intensity per count in RGB color space")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file name")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (outf.empty() || tsvf.empty() || cmapf.empty() || minmaxf.empty())
        error("--tsv, --minmax, --cmap and --out must be specified");

    notice("Analysis started");

    // read the manifest file and determine the xmin/xmax/ymin/ymax
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    read_minmax_double(minmaxf.c_str(), xmin, xmax, ymin, ymax);
    notice("Bounding box: xmin=%llu, xmax=%llu, ymin=%llu, ymax=%llu", 
        (unsigned long long)xmin, (unsigned long long)xmax, 
        (unsigned long long)ymin, (unsigned long long)ymax);
 
    // parse the color map file
    tsv_reader cmap_tr(cmapf.c_str());
    int32_t icol_gene = -1, icol_hex = -1;
    std::map<std::string, uint32_t> gene2rgb;
    std::map<std::string, uint32_t>::iterator gene2rgb_it;
    while (cmap_tr.read_line() > 0)
    {
        if ( icol_gene < 0 ) {
            for (int32_t i = 0; i < cmap_tr.nfields; ++i)
            {
                if (strcmp(cmap_tr.str_field_at(i), cmap_colname_gene.c_str()) == 0)
                    icol_gene = i;
                else if (strcmp(cmap_tr.str_field_at(i), cmap_colname_hex.c_str()) == 0)
                    icol_hex = i;
            }
            if (icol_gene < 0 || icol_hex < 0)
                error("Cannot find column %s or %s in the color map file %s",
                      cmap_colname_gene.c_str(), cmap_colname_hex.c_str(), cmapf.c_str());
        }
        else {
            const char* gene = cmap_tr.str_field_at(icol_gene);
            const char* hex = cmap_tr.str_field_at(icol_hex);
            uint32_t rgba = hex2rgba(hex);
            if ( gene2rgb.find(gene) != gene2rgb.end() )
            {
                error("Gene %s is already defined in the color map.");
            }
            gene2rgb[gene] = rgba;
        }
    }
    cmap_tr.close();

    int32_t height = (int32_t)(ceil((double)(ymax - ymin + 1) / coord_per_pixel));
    int32_t width = (int32_t)(ceil((double)(xmax - xmin + 1) / coord_per_pixel));
    cimg_library::CImg<unsigned char> image(width, height, 1, 3, 0);

    // read the SGE matrix
    notice("Processing TSV file %s ...", tsvf.c_str());
    tsv_reader tsv_tr(tsvf.c_str());
    int32_t icol_x = -1, icol_y = -1, icol_count = -1;
    icol_gene = -1;
    if ( skip_tsv_header ) {
        icol_x = 0; icol_y = 1; icol_gene = 2; icol_count = 3;
    }
    uint64_t nlines = 0, npass = 0, nskip = 0;
    while ( tsv_tr.read_line() ) 
    {
        if ( icol_x < 0 ) {
            for (int32_t i = 0; i < tsv_tr.nfields; ++i)
            {
                if (strcmp(tsv_tr.str_field_at(i), tsv_colname_X.c_str()) == 0)
                    icol_x = i;
                else if (strcmp(tsv_tr.str_field_at(i), tsv_colname_Y.c_str()) == 0)
                    icol_y = i;
                else if (strcmp(tsv_tr.str_field_at(i), tsv_colname_gene.c_str()) == 0)
                    icol_gene = i;
                else if (strcmp(tsv_tr.str_field_at(i), tsv_colname_count.c_str()) == 0)
                    icol_count = i;
            }
            if (icol_x < 0 || icol_y < 0 || icol_gene < 0 || icol_count < 0)
                error("Cannot find column %s, %s, %s or %s in the TSV file %s",
                      tsv_colname_X.c_str(), tsv_colname_Y.c_str(),
                      tsv_colname_gene.c_str(), tsv_colname_count.c_str(), tsvf.c_str());
        }
        else {
            const char* gene = tsv_tr.str_field_at(icol_gene);
            gene2rgb_it = gene2rgb.find(gene);
            if (gene2rgb_it == gene2rgb.end())
            {
                ++nskip;
            }
            else {
                double x = tsv_tr.double_field_at(icol_x);
                double y = tsv_tr.double_field_at(icol_y);
                int32_t px = (int32_t)floor((x - xmin) / coord_per_pixel);
                int32_t py = (int32_t)floor((y - ymin) / coord_per_pixel);
                for(int32_t j=0; j < 3; ++j) {
                    int32_t c = image(px, py, j);
                    int32_t rgb = (gene2rgb_it->second >> (8 * (2 - j))) & 0xFF;
                    c += (int32_t)(rgb * intensity_per_count * tsv_tr.double_field_at(icol_count));
                    if (c > 255) c = 255;
                    image(px, py, j) = c;
                }
                ++npass;
            }
            
            ++nlines;
            if (nlines % 10000000 == 0)
                notice("Processing %llu lines (%llu pass, %llu skip) from %s", nlines, npass, nskip, tsvf.c_str());
        }
    }
    tsv_tr.close();
    notice("Finished processin g %llu lines (%llu pass, %llu skip) from %s", nlines, npass, nskip, tsvf.c_str());

    notice("Writing the image to %s", outf.c_str());
    image.save_png(outf.c_str());

    return 0;
}
