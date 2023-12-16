#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include <cmath>
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#define cimg_display 0  // remove the need for X11 library
#include "cimg/CImg.h"

/////////////////////////////////////////////////////////////////////////
// draw-xy : Draw the single-color image of points in 2D space
////////////////////////////////////////////////////////////////////////
int32_t cmdDraw3way(int32_t argc, char **argv)
{
    std::string manifestf;
    std::string nbcdf;
    std::vector<std::string> nmatchfs;
    std::string ngebcdf;
    double coord_per_pixel = 1000.; // 1 pixel = 1000 units = 1 um
    int32_t icol_x_nbcd = 3;
    int32_t icol_y_nbcd = 4;
    int32_t icol_x_nmatch = 3;
    int32_t icol_y_nmatch = 4;
    int32_t icol_x_ngebcd = 5;
    int32_t icol_y_ngebcd = 6;
    int32_t icol_gene_ngebcd = 7;
    int32_t isubcol_gene_ngebcd = 0;

    std::string color_nbcd = "#320000";
    std::string color_nmatch = "#000032";
    std::string color_ngebcd = "#006400";
    std::string id_manifest = "1_1";

    std::string outf;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input files", NULL)
    LONG_STRING_PARAM("manifest", &manifestf, "Manifest file from combine-sbcds. xmin/xmax/ymin/ymax will be automatically detected")
    LONG_STRING_PARAM("nbcd", &nbcdf, "Spatial barcode dictionary generated from 'combine-sbcds' command")
    LONG_MULTI_STRING_PARAM("nmatch", &nmatchfs, "Spatial barcode dictionary generated from 'match-sbcds' command")
    LONG_STRING_PARAM("ngebcd", &ngebcdf, "Spatial barcode dictionary generated from alignment pipeline")

    LONG_PARAM_GROUP("Input options", NULL)
    LONG_INT_PARAM("icol-x-nbcd", &icol_x_nbcd, "0-based index of the column for x in nbcd")
    LONG_INT_PARAM("icol-y-nbcd", &icol_y_nbcd, "0-based index of the column for y in nbcd")
    LONG_INT_PARAM("icol-x-nmatch", &icol_x_nmatch, "0-based index of the column for x in nmatch")
    LONG_INT_PARAM("icol-y-nmatch", &icol_y_nmatch, "0-based index of the column for y in nmatch")
    LONG_INT_PARAM("icol-x-ngebcd", &icol_x_ngebcd, "0-based index of the column for x in ngebcd")
    LONG_INT_PARAM("icol-y-ngebcd", &icol_y_ngebcd, "0-based index of the column for y in ngebcd")
    LONG_INT_PARAM("icol-gene-ngebcd", &icol_gene_ngebcd, "0-based index of the column for gene count in ngebcd")
    LONG_INT_PARAM("isubcol-gene-ngebcd", &isubcol_gene_ngebcd, "0-based index of the sub column for gene count in ngebcd (0: Gene, 1: GeneFull, ...)")
    LONG_STRING_PARAM("id-manifest", &id_manifest, "ID of the tile in the manifest file to use")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_DOUBLE_PARAM("coord-per-pixel", &coord_per_pixel, "Number of coordinate units per pixel")
    LONG_STRING_PARAM("color-nbcd", &color_nbcd, "RGB hex color code for nbcd per observation")
    LONG_STRING_PARAM("color-nmatch", &color_nbcd, "RGB hex color code for nmatch per observation")
    LONG_STRING_PARAM("color-nge", &color_nbcd, "RGB hex color code for nge per observation")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file name")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( manifestf.empty() || outf.empty() )
        error("--manifest and --out must be specified");

    notice("Analysis started");

    // read the manifest file and determine the xmin/xmax/ymin/ymax
    uint64_t xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    dataframe_t df(manifestf.c_str());
    if ( df.nrows == 0 )
        error("Empty dataframe %s", manifestf.c_str());
    int32_t i_id = df.get_colidx("id");
    if ( i_id < 0 )
        error("Cannot find the column 'id' in %s", manifestf.c_str());
    for(int32_t i=0; df.nrows; ++i) {
        // find the matching tile
        if ( df.get_str_elem(i, i_id).compare(id_manifest) == 0 ) {
            // parse xmin/xmax/ymin/ymax
            int32_t i_xmin = df.get_colidx("xmin");
            int32_t i_xmax = df.get_colidx("xmax");
            int32_t i_ymin = df.get_colidx("ymin");
            int32_t i_ymax = df.get_colidx("ymax");
            if ( i_xmin < 0 || i_xmax < 0 || i_ymin < 0 || i_ymax < 0 )
                error("Cannot find the columns 'xmin', 'xmax', 'ymin', 'ymax' in %s", manifestf.c_str());
            xmin = df.get_uint64_elem(i, i_xmin);
            xmax = df.get_uint64_elem(i, i_xmax);
            ymin = df.get_uint64_elem(i, i_ymin);
            ymax = df.get_uint64_elem(i, i_ymax);
            break;
        }
    }

    // parse input colors
    uint32_t rgb_nbcd[3], rgb_nmatch[3], rgb_ngebcd[3];
    sscanf(color_nbcd.c_str() + 1, "%02x%02x%02x", &rgb_nbcd[0], &rgb_nbcd[1], &rgb_nbcd[2]);
    sscanf(color_nmatch.c_str() + 1, "%02x%02x%02x", &rgb_nmatch[0], &rgb_nmatch[1], &rgb_nmatch[2]);
    sscanf(color_ngebcd.c_str() + 1, "%02x%02x%02x", &rgb_ngebcd[0], &rgb_ngebcd[1], &rgb_ngebcd[2]);

    uint64_t width  = (uint64_t)(ceil((double)(xmax - xmin + 1.0) / (double)coord_per_pixel));
    uint64_t height = (uint64_t)(ceil((double)(ymax - ymin + 1.0) / (double)coord_per_pixel));
    if ( width == 0 || height == 0 )
        error("Invalid width/height: %llu/%llu", width, height);

    // create and image
    cimg_library::CImg<unsigned char> image(width, height, 1, 3, 255);

    // read the nbcd file
    if ( ! nbcdf.empty() ) {
        tsv_reader nbcd(nbcdf.c_str());
        while ( nbcd.read_line() ) {
            if ( nbcd.nfields <= icol_x_nbcd || nbcd.nfields <= icol_y_nbcd )
                error("Input file %s does not have enough columns - only %d", nbcdf.c_str(), nbcd.nfields);

            uint64_t x = nbcd.uint64_field_at(icol_x_nbcd);
            uint64_t y = nbcd.uint64_field_at(icol_y_nbcd);
            int32_t ix = (int32_t)((x - xmin)/ coord_per_pixel);
            int32_t iy = (int32_t)((y - ymin) / coord_per_pixel);
            if ( ix > width || iy > height )
                error("Out of range point detected (%lf, %lf): width = %llu, height = %llu, coord_per_pixel = %lf", x, y, width, height, coord_per_pixel);

            for(int32_t i=0; i < 3; ++i) {
                int32_t c = image(ix, iy, i);
                if ( c + rgb_nbcd[i] > 255 ) c = 255;
                else c += rgb_nbcd[i];
                image(ix, iy, i) = c;
            }
        }
    }

    // read the nmatch files
    for(int32_t i=0; i < nmatchfs.size(); ++i) {
        tsv_reader tf(nmatchfs[i].c_str());
        while ( tf.read_line() ) {
            if ( tf.nfields <= icol_x_nmatch || tf.nfields <= icol_y_nmatch )
                error("Input file %s does not have enough columns - only %d", nmatchfs[i].c_str(), tf.nfields);

            uint64_t x = tf.uint64_field_at(icol_x_nmatch);
            uint64_t y = tf.uint64_field_at(icol_y_nmatch);
            int32_t ix = (int32_t)((x - xmin)/ coord_per_pixel);
            int32_t iy = (int32_t)((y - ymin) / coord_per_pixel);
            if ( ix > width || iy > height )
                error("Out of range point detected (%lf, %lf): width = %llu, height = %llu, coord_per_pixel = %lf", x, y, width, height, coord_per_pixel);

            for(int32_t i=0; i < 3; ++i) {
                int32_t c = image(ix, iy, i);
                if ( c + rgb_nmatch[i] > 255 ) c = 255;
                else c += rgb_nmatch[i];
                image(ix, iy, i) = c;
            }
        }
    }

    // read the nge barcode files
    if ( ! ngebcdf.empty() ) {
        tsv_reader tf(ngebcdf.c_str());
        while ( tf.read_line() ) {
            if ( tf.nfields <= icol_x_ngebcd || tf.nfields <= icol_y_ngebcd || tf.nfields <= icol_gene_ngebcd )
                error("Input file %s does not have enough columns - only %d", ngebcdf.c_str(), tf.nfields);

            double x = tf.double_field_at(icol_x_ngebcd);
            double y = tf.double_field_at(icol_y_ngebcd);        
            int32_t ix = (int32_t)((x - xmin)/ coord_per_pixel);
            int32_t iy = (int32_t)((y - ymin) / coord_per_pixel);
            const char* s = tf.str_field_at(icol_gene_ngebcd);
            for(int32_t i=0; i < isubcol_gene_ngebcd; ++i) {
                s = strchr(s, ',') + 1;
            }
            int32_t gene_count = atoi(s);

            if ( ix > width || iy > height )
                error("Out of range point detected (%lf, %lf): width = %llu, height = %llu, coord_per_pixel = %lf", x, y, width, height, coord_per_pixel);

            for(int32_t i=0; i < 3; ++i) {
                int32_t c = image(ix, iy, i);
                if ( c + rgb_ngebcd[i] * gene_count > 255 ) c = 255;
                else c += ( rgb_ngebcd[i] * gene_count );
                image(ix, iy, i) = c;
            }
        }
    }

    image.save_png(outf.c_str());

    return 0;
}
