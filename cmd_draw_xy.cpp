#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "sge.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
//#include <gd.h>
#define cimg_display 0  // remove the need for X11 library
#include "cimg/CImg.h"

/////////////////////////////////////////////////////////////////////////
// draw-xy : Draw the image of points in 2D space
////////////////////////////////////////////////////////////////////////
int32_t cmdDrawXY(int32_t argc, char **argv)
{
    std::string tsvf;
    int32_t icolx = 0;
    int32_t icoly = 1;
    double coord_per_pixel = 1.0; // 1 pixel = 1 coordinate unit
    int32_t width = 0;
    int32_t height = 0;
    int32_t intensity_per_obs = 1; // intensity per observation
    std::string outf;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("tsv", &tsvf, "tsv file to draw the x-y coordinates. /dev/stdin for stdin")
    LONG_INT_PARAM("icol-x", &icolx, "0-based index of the column for x")
    LONG_INT_PARAM("icol-y", &icoly, "0-based index of the column for y")

    LONG_PARAM_GROUP("Settings", NULL)
    LONG_INT_PARAM("width", &width, "Width of the image")
    LONG_INT_PARAM("height", &height, "Height of the image")
    LONG_DOUBLE_PARAM("coord-per-pixel", &coord_per_pixel, "Number of coordinate units per pixel")
    LONG_INT_PARAM("intensity-per-obs", &intensity_per_obs, "Intensity per pixel per observation")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file name")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( tsvf.empty() || outf.empty() )
        error("--tsv and --out must be specified");
    if ( width == 0 || height == 0 )
        error("--width and --height must be specified");

    notice("Analysis started");

    //gdImagePtr im = gdImageCreateTrueColor(width, height);
    cimg_library::CImg<unsigned char> image(width, height, 1, 3, 255);

    //int32_t colors[256];
    //for(int32_t i=0; i < 256; ++i) 
    //    colors[i] = gdImageColorAllocate(im, i, i, i);

    uint8_t *imbuf = (uint8_t *)malloc(width * height);
    memset(imbuf, 0, width * height);

    tsv_reader tf(tsvf.c_str());

//    FILE *pngout = fopen(outf.c_str(), "wb");
//    if (!pngout) {
//        error("Cannot open %s for writing", outf.c_str());
 //   }

    while ( tf.read_line() ) {
        if ( tf.nfields <= icolx || tf.nfields <= icoly )
            error("Input file %s does not have enough columns - only %d", tsvf.c_str(), tf.nfields);

        double x = tf.double_field_at(icolx);
        double y = tf.double_field_at(icoly);
        int32_t ix = (int32_t)(x / coord_per_pixel);
        int32_t iy = (int32_t)(y / coord_per_pixel);
        if ( ix >= 0 && ix < width && iy >= 0 && iy < height ) {
            if ( imbuf[iy * width + ix] + intensity_per_obs < 256 )
                imbuf[iy * width + ix] += (uint8_t)intensity_per_obs;
            else
                imbuf[iy * width + ix] = 255;
        }
        else {
            error("Out of range point detected (%lf, %lf)", x, y);
        }
    }

    for(int32_t iy=0; iy < height; ++iy) {
        for(int32_t ix=0; ix < width; ++ix) {
            //gdImageSetPixel(im, ix, iy, colors[imbuf[iy * width + ix]]);
            uint8_t c = imbuf[iy * width + ix];
            image(ix, iy, 0) = c;
            image(ix, iy, 1) = c;
            image(ix, iy, 2) = c;
        }
    }

    image.save_png(outf.c_str());
    free(imbuf);

    return 0;
}
