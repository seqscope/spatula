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
// draw-xy : Draw the single-color image of points in 2D space
////////////////////////////////////////////////////////////////////////
int32_t cmdDrawXY(int32_t argc, char **argv)
{
    std::string tsvf;
    int32_t icolx = 0;
    int32_t icoly = 1;
    double coord_per_pixel = 1.0; // 1 pixel = 1 coordinate unit
    int32_t width = 0;
    int32_t height = 0;
    int32_t intensity_per_obs = 1;  // intensity per observation
    int32_t verbose_freq = 1000000; // report frequency of input reading
    bool auto_adjust_intensity = false;
    int32_t max_intensity = 255;
    double auto_adjust_quantile = 0.99;
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
    LONG_PARAM("auto-adjust", &auto_adjust_intensity, "Automatically adjust the intensity of the color based on the maximum count")
    LONG_DOUBLE_PARAM("adjust-quantile", &auto_adjust_quantile, "Quantile of pixel to use for auto-adjustment among non-zero pixels")
    LONG_INT_PARAM("max-intensity", &max_intensity, "Maximum value of possible intensity")


    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output file name")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( tsvf.empty() || outf.empty() )
        error("--tsv and --out must be specified");
    if ( width == 0 || height == 0 ) 
        notice("width and height is not specified. Will be automatically detected while reading the input files");

    notice("Analysis started");

    std::vector<uint8_t*> imbufs;
    int32_t cur_height = height == 0 ? 1000 : height;
    int32_t max_y = 0;

    tsv_reader tf(tsvf.c_str());

    std::vector<int32_t> intensity_counts(max_intensity + 1, 0);

    uint64_t nlines = 0;
    while ( tf.read_line() ) {
        if ( tf.nfields <= icolx || tf.nfields <= icoly )
            error("Input file %s does not have enough columns - only %d", tsvf.c_str(), tf.nfields);

        double x = tf.double_field_at(icolx);
        double y = tf.double_field_at(icoly);
        int32_t ix = (int32_t)(x / coord_per_pixel);
        int32_t iy = (int32_t)(y / coord_per_pixel);

        // if width and height are specified, check if the point is within the range
        if ( ( width > 0 ) && ( ix < 0 || ix >= width ) ) {
            error("Out of range point detected (%lf, %lf): width = %d, coord_per_pixel = %lf", x, y, width, coord_per_pixel);
        }
        if ( ( height > 0 ) && ( iy < 0 || iy >= height ) ) {
            error("Out of range point detected (%lf, %lf): width = %d, coord_per_pixel = %lf", x, y, width, coord_per_pixel);
        }

        // if the x-coordinate is out of range, add more coordinates
        if ( ix >= imbufs.size() ) imbufs.resize(ix + 1, NULL);
        if ( imbufs[ix] == NULL ) {
            imbufs[ix] = (uint8_t *)calloc(cur_height, sizeof(uint8_t));
        }

        // if the y-coordinate is out of range, double the cur_height
        if ( cur_height <= iy ) {
            int32_t new_height = cur_height;
            while ( new_height <= iy ) new_height *= 2;
            for(int32_t i=0; i < (int32_t)imbufs.size(); ++i) {
                imbufs[i] = (uint8_t *)realloc(imbufs[i], new_height * sizeof(uint8_t));
                memset(imbufs[i] + cur_height, 0, (new_height - cur_height) * sizeof(uint8_t));
            }
            cur_height = new_height;
        }

        max_y = iy > max_y ? iy : max_y;

        if ( imbufs[ix][iy] + intensity_per_obs <= max_intensity ) {
            if ( imbufs[ix][iy] == 0 ) {
                intensity_counts[intensity_per_obs]++;
            }
            else {
                intensity_counts[imbufs[ix][iy]]--;
                intensity_counts[imbufs[ix][iy] + intensity_per_obs]++;
            }
            imbufs[ix][iy] += (uint8_t)intensity_per_obs;
        }
        else {
            if ( imbufs[ix][iy] == 0 ) {
                intensity_counts[intensity_per_obs]++;
            }
            else {
                intensity_counts[imbufs[ix][iy]]--;
                intensity_counts[max_intensity]++;
            }
            imbufs[ix][iy] = max_intensity;    
        }

	++nlines;

	if ( nlines % verbose_freq == 0 ) {
	    notice("Reading %llu input lines... max_y = %d, cur_height = %d", nlines, max_y, cur_height);
	}
    }

    if ( width == 0 ) {
        width = (int32_t)imbufs.size();
        notice("Setting the width = %d", width);
    }
    if ( height == 0 ) {
        height = max_y + 1;
        notice("Setting the height = %d", height);
    }

    if ( auto_adjust_intensity ) {
        notice("Auto-adjusting the intensity of the image");

        // obtain quantile threshold of non-zero intensities
        std::vector<int32_t> cumulative_counts(max_intensity + 1, 0);
        for(int32_t i=1; i <= max_intensity; ++i) {
            cumulative_counts[i] = cumulative_counts[i-1] + intensity_counts[i];
        }

        int32_t threshold_intensity = max_intensity;
        for(int32_t i=max_intensity; i > 0; --i) {
            if ( cumulative_counts[i] > auto_adjust_quantile * cumulative_counts[max_intensity] ) {
                threshold_intensity = i;
            }
            else {
                break;
            }
        }

        notice("Auto-adjusted intensity threshold = %d", threshold_intensity);
        for(int32_t ix=0; ix < width; ++ix) {
            if ( imbufs[ix] != NULL ) {
                for(int32_t iy=0; iy < height; ++iy) {
                    if ( imbufs[ix][iy] > 0 ) {
                        if ( imbufs[ix][iy] > threshold_intensity ) {
                            imbufs[ix][iy] = max_intensity;
                        }
                        else {
                            imbufs[ix][iy] = (int32_t)(imbufs[ix][iy] * max_intensity / threshold_intensity);
                        }
                    }
                }
            }
        }
    }

    notice("Creating an image in memory");
    cimg_library::CImg<unsigned char> image(width, height, 1, 3, 0);
    
    for(int32_t ix=0; ix < width; ++ix) {
        if ( imbufs[ix] == NULL ) {
            for(int32_t iy=0; iy < height; ++iy) {
                image(ix, iy, 0) = 0;
                image(ix, iy, 1) = 0;
                image(ix, iy, 2) = 0;
            }
        }
        else {
            for(int32_t iy=0; iy < height; ++iy) {
                uint8_t c = imbufs[ix][iy];
                image(ix, iy, 0) = c;
                image(ix, iy, 1) = c;
                image(ix, iy, 2) = c;
            }
        }
    }

    notice("Saving the png file to %s", outf.c_str());
    image.save_png(outf.c_str());

    notice("Freeing up the memory..");
    for(int32_t ix=0; ix < width; ++ix) {
        if ( imbufs[ix] != NULL ) 
            free(imbufs[ix]);
    }

    return 0;
}
