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
#include <cmath>
#define cimg_display 0  // remove the need for X11 library
#include "cimg/CImg.h"
//#include <png.h>

// Define FPNG_USE_ZLIB 1 before including fpng.h to enable zlib for robust decoding.
// This is highly recommended if your input PNGs can come from various sources.
#define FPNG_USE_ZLIB 1
#include "fpng/fpng.h"

// Function to load a grayscale PNG file
// Returns true on success, false on failure.
// Fills width, height, and populates grayscale_data (vector of row pointers).
// grayscale_data will contain raw pixel values (one byte per pixel).
// bool load_grayscale_png(const char* filename, int& width, int& height, std::vector<png_bytep>& row_pointers_gray, png_byte& bit_depth_param) {
//     FILE *fp = fopen(filename, "rb");
//     if (!fp) {
//         error("Could not open file %s for reading.", filename);
//     }

//     png_structp png_ptr_read = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//     if (!png_ptr_read) {
//         error("png_create_read_struct failed.");
//         fclose(fp);
//         return false;
//     }

//     png_infop info_ptr_read = png_create_info_struct(png_ptr_read);
//     if (!info_ptr_read) {
//         error("png_create_info_struct failed.");
//         png_destroy_read_struct(&png_ptr_read, (png_infopp)NULL, (png_infopp)NULL);
//         fclose(fp);
//         return false;
//     }

//     if (setjmp(png_jmpbuf(png_ptr_read))) {
//         error("Error during png_init_io or png_read_info.");
//         png_destroy_read_struct(&png_ptr_read, &info_ptr_read, (png_infopp)NULL);
//         fclose(fp);
//         return false;
//     }

//     // Initialize PNG I/O
//     png_init_io(png_ptr_read, fp);
//     png_read_info(png_ptr_read, info_ptr_read);

//     width = png_get_image_width(png_ptr_read, info_ptr_read);
//     height = png_get_image_height(png_ptr_read, info_ptr_read);
//     png_byte original_color_type = png_get_color_type(png_ptr_read, info_ptr_read);
//     png_byte original_bit_depth = png_get_bit_depth(png_ptr_read, info_ptr_read);

//     // 1. Convert 16-bit per channel to 8-bit per channel
//     if (original_bit_depth == 16) {
//         png_set_strip_16(png_ptr_read);
//     }

//     // 2. Convert palette images to RGB
//     if (original_color_type == PNG_COLOR_TYPE_PALETTE) {
//         png_set_palette_to_rgb(png_ptr_read);
//         // Image is now effectively RGB for further processing
//         original_color_type = PNG_COLOR_TYPE_RGB; // Update effective type
//         original_bit_depth = 8; // Palette to RGB results in 8-bit per channel
//     }

//     // 3. Convert RGB or RGBA to Grayscale
//     // This handles your specific issue with Color Type 2 (RGB)
//     if (original_color_type == PNG_COLOR_TYPE_RGB || original_color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
//         // Parameters for rgb_to_gray: error_action, red_weight, green_weight
//         // Using -1 for weights uses libpng's default calculation (ITU-R BT.709 Luma)
//         png_set_rgb_to_gray(png_ptr_read, 1, -1, -1);
//         // Update effective type: if it was RGBA, it's now GRAY_ALPHA, else GRAY
//         original_color_type = (original_color_type == PNG_COLOR_TYPE_RGB_ALPHA) ?
//                               PNG_COLOR_TYPE_GRAY_ALPHA : PNG_COLOR_TYPE_GRAY;
//     }

//     // 4. If image is Grayscale with Alpha, strip the alpha channel to get pure grayscale
//     if (original_color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
//         png_set_strip_alpha(png_ptr_read);
//         original_color_type = PNG_COLOR_TYPE_GRAY; // Now it's pure GRAY
//     }

//     // 5. Expand low bit-depth grayscale (1, 2, 4 bits) to 8 bits per pixel
//     // This should be applied if the image is now effectively GRAY and had an original low bit depth
//     // Note: png_get_bit_depth before png_read_update_info still gives original depth.
//     // We only do this if the image is (now) GRAY and its original bit depth was < 8.
//     png_byte current_bit_depth = png_get_bit_depth(png_ptr_read, info_ptr_read); // Gets original bit depth
//     if (original_color_type == PNG_COLOR_TYPE_GRAY && current_bit_depth < 8) {
//         png_set_expand_gray_1_2_4_to_8(png_ptr_read); // Preferred modern way
//         // png_set_packing(png_ptr_read); // Older way, also works
//     }

//     // Apply all requested transformations
//     png_read_update_info(png_ptr_read, info_ptr_read);

//     // ---- Validation after transformations ----
//     // The image data should now be 8-bit grayscale
//     png_byte final_color_type = png_get_color_type(png_ptr_read, info_ptr_read);
//     png_byte final_bit_depth = png_get_bit_depth(png_ptr_read, info_ptr_read);

//     bit_depth_param = final_bit_depth; // Update the output parameter

//     if (final_color_type != PNG_COLOR_TYPE_GRAY || final_bit_depth != 8) {
//         notice("Original color type: %d, Original bit depth: %d", original_color_type, original_bit_depth);
//         notice("After transformations - Color type: %d (expected %d)", final_color_type, PNG_COLOR_TYPE_GRAY);
//         notice("After transformations - Bit depth: %d (expected 8)", final_bit_depth);
//         error("Failed to convert input PNG to 8-bit grayscale.");
//         png_destroy_read_struct(&png_ptr_read, &info_ptr_read, (png_infopp)NULL);
//         fclose(fp);
//         return false;
//     }

//     if (setjmp(png_jmpbuf(png_ptr_read))) {
//         error("Error during png_read_image.");
//         png_destroy_read_struct(&png_ptr_read, &info_ptr_read, (png_infopp)NULL);
//         if (!row_pointers_gray.empty()) {
//             for (int y = 0; y < height; y++) {
//                 if(row_pointers_gray[y]) free(row_pointers_gray[y]);
//             }
//         }
//         fclose(fp);
//         return false;
//     }

//     row_pointers_gray.resize(height);
//     for (int y = 0; y < height; y++) {
//         row_pointers_gray[y] = (png_bytep)malloc(png_get_rowbytes(png_ptr_read, info_ptr_read));
//         if (!row_pointers_gray[y]) {
//             error("Could not allocate memory for grayscale row %d", y);
//             for(int k=0; k<y; ++k) free(row_pointers_gray[k]);
//             row_pointers_gray.clear();
//             png_destroy_read_struct(&png_ptr_read, &info_ptr_read, (png_infopp)NULL);
//             fclose(fp);
//             return false;
//         }
//     }

//     png_read_image(png_ptr_read, row_pointers_gray.data());

//     fclose(fp);
//     png_destroy_read_struct(&png_ptr_read, &info_ptr_read, (png_infopp)NULL);
//     notice("Grayscale PNG loaded successfully: %dx%d @ %dbpp", width, height, (int)bit_depth_param);
//     return true;
// }

// // Function to save RGBA data as a PNG file
// // rgba_data should be a flat vector: [R,G,B,A, R,G,B,A, ...]
// // Returns true on success, false on failure.
// bool save_rgba_png(const char* filename, int width, int height, const std::vector<unsigned char>& rgba_data) {
//     if (rgba_data.size() != static_cast<size_t>(width * height * 4)) {
//         error("Error: RGBA data size does not match image dimensions (%d x %d)", width, height);
//         return false;
//     }

//     FILE *fp = fopen(filename, "wb");
//     if (!fp) {
//         error("Error: Could not open file %s for writing.", filename);
//         return false;
//     }

//     png_structp png_ptr_write = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//     if (!png_ptr_write) {
//         error("Error: png_create_write_struct failed.");
//         fclose(fp);
//         return false;
//     }

//     png_infop info_ptr_write = png_create_info_struct(png_ptr_write);
//     if (!info_ptr_write) {
//         error("Error: png_create_info_struct failed.");
//         png_destroy_write_struct(&png_ptr_write, (png_infopp)NULL);
//         fclose(fp);
//         return false;
//     }

//     if (setjmp(png_jmpbuf(png_ptr_write))) {
//         error("Error: Error during png_init_io or png_set_compression_level.");
//         png_destroy_write_struct(&png_ptr_write, &info_ptr_write);
//         fclose(fp);
//         return false;
//     }

//     png_init_io(png_ptr_write, fp);

//     // Set PNG header information
//     png_set_IHDR(
//         png_ptr_write,
//         info_ptr_write,
//         width,
//         height,
//         8, // Bit depth per channel (RGBA usually 8)
//         PNG_COLOR_TYPE_RGBA,
//         PNG_INTERLACE_NONE,
//         PNG_COMPRESSION_TYPE_DEFAULT,
//         PNG_FILTER_TYPE_DEFAULT
//     );
//     png_write_info(png_ptr_write, info_ptr_write);

//     // Create row pointers for the RGBA data
//     std::vector<png_bytep> row_pointers_rgba(height);
//     for (int y = 0; y < height; y++) {
//         // const_cast is safe here as libpng won't modify input data during write
//         row_pointers_rgba[y] = const_cast<png_bytep>(rgba_data.data() + y * width * 4);
//     }

//     if (setjmp(png_jmpbuf(png_ptr_write))) {
//         error("Error: Error during png_write_image.");
//         png_destroy_write_struct(&png_ptr_write, &info_ptr_write);
//         fclose(fp);
//         return false;
//     }

//     png_write_image(png_ptr_write, row_pointers_rgba.data());

//     if (setjmp(png_jmpbuf(png_ptr_write))) {
//         error("Error: Error during png_write_end.");
//         png_destroy_write_struct(&png_ptr_write, &info_ptr_write);
//         fclose(fp);
//         return false;
//     }
//     png_write_end(png_ptr_write, NULL);

//     fclose(fp);
//     png_destroy_write_struct(&png_ptr_write, &info_ptr_write);
//     notice("RGBA PNG saved successfully: %s", filename);
//     return true;
// }

// /////////////////////////////////////////////////////////////////////////
// // png-mono2rgba : Convert a monochrome PNG image to RGBA format
// ////////////////////////////////////////////////////////////////////////
// int32_t cmdPngMono2RGBA(int32_t argc, char **argv)
// {
//     std::string in_png;  // input PNG file (grayscale)
//     std::string out_png; // output PNG file
//     int32_t thres_transparency = 0;  // threshold for transparent pixels
//     std::string rgb("FFFFFF");  // RGB color for the monochrome image

//     paramList pl;
//     BEGIN_LONG_PARAMS(longParameters)
//     LONG_PARAM_GROUP("Input options", NULL)
//     LONG_STRING_PARAM("in", &in_png, "Input PNG file (grayscale)")
//     LONG_STRING_PARAM("out", &out_png, "Output PNG file")
//     LONG_INT_PARAM("thres", &thres_transparency, "Threshold for transparent pixels (default: 0)")
//     LONG_STRING_PARAM("rgb", &rgb, "RGB color for the monochrome image (default: #FFFFFF)")
//     END_LONG_PARAMS();

//     pl.Add(new longParams("Available Options", longParameters));
//     pl.Read(argc, argv);
//     pl.Status();

//     notice("Analysis started");
//     if ( in_png.empty() || out_png.empty() )
//         error("--in and --out must be specified");

//     const char* rgb_c = rgb.c_str();
//     if ( rgb_c[0] == '#' ) ++rgb_c; // remove the # sign
//     uint32_t rgb_i[3];
//     sscanf(rgb_c, "%02x%02x%02x", &rgb_i[0], &rgb_i[1], &rgb_i[2]); // parse RGB hex code
    

//     int32_t width, height;
//     png_byte bit_depth_gray;
//     std::vector<png_bytep> grayscale_rows; // Vector of row pointers

//     notice("Loading %s", in_png.c_str());
//     if (!load_grayscale_png(in_png.c_str(), width, height, grayscale_rows, bit_depth_gray)) {
//         return 1;
//     }

//     // Convert grayscale to RGBA
//     std::vector<unsigned char> rgba_data;
//     rgba_data.reserve(width * height * 4);

//     notice("Processing individual pixels");
//     for (int32_t y = 0; y < height; y++) {
//         png_bytep row = grayscale_rows[y];
//         for (int32_t x = 0; x < width; x++) {
//             png_byte gray_value = row[x]; // Assuming 1 byte per pixel for grayscale
//             if (gray_value <= thres_transparency) {
//                 rgba_data.push_back(0); // R (transparent)
//                 rgba_data.push_back(0); // G
//                 rgba_data.push_back(0); // B
//                 rgba_data.push_back(0); // Alpha
//             }
//             else {            
//                 rgba_data.push_back((png_byte)floor(gray_value / 255.0 * rgb_i[0])); // R
//                 rgba_data.push_back((png_byte)floor(gray_value / 255.0 * rgb_i[1])); // R
//                 rgba_data.push_back((png_byte)floor(gray_value / 255.0 * rgb_i[2])); // R
//                 rgba_data.push_back(gray_value); // Alpha
//             }
//         }
//     }

//     notice("Freeing unnecessary memory");
//     for (int y = 0; y < height; y++) {
//         if (grayscale_rows[y]) {
//             free(grayscale_rows[y]);
//         }
//     }
//     grayscale_rows.clear();


//     notice("Saving the png file to %s", out_png.c_str());
//     if (!save_rgba_png(out_png.c_str(), width, height, rgba_data)) {
//         return 1;
//     }
//     return 0;
// }


// Function to save RGBA data as a PNG file using FPNG
bool save_rgba_with_fpng(const char* filename,
                         const std::vector<uint8_t>& rgba_data,
                         uint32_t width,
                         uint32_t height) {
    if (rgba_data.empty() || rgba_data.size() != static_cast<size_t>(width * height * 4)) {
        error("Error: RGBA data is empty or has incorrect size for saving.");
        return false;
    }

    // FPNG_ENCODE_SLOWER_BUT_SMALLER can be used for better compression, or 0 for defaults.
    if (!fpng::fpng_encode_image_to_file(filename, rgba_data.data(), width, height, 4 /* num_chans for RGBA */)) {
        error("Error: FPNG failed to encode image to file: %s", filename);
        // FPNG encoding doesn't provide a status code directly in this function signature,
        // but failure is indicated by the boolean return.
        return false;
    }

    return true;
}


/////////////////////////////////////////////////////////////////////////
// png-mono2rgba : Convert a monochrome PNG image to RGBA format
////////////////////////////////////////////////////////////////////////
int32_t cmdPngMono2RGBA(int32_t argc, char **argv)
{
    std::string in_png;  // input PNG file (grayscale)
    std::string out_png; // output PNG file
    int32_t thres_transparency_below = 0;  // threshold for transparent pixels
    int32_t thres_transparency_above = 255;       // threshold for transparent pixels
    std::string rgb("FFFFFF");  // RGB color for the monochrome image
    double alpha_linear_coef = 1.0;
    int32_t alpha_constant_max = 255;
    double alpha_quadratic_coef = 16.0;   
    bool use_quadratic_alpha = false; // Use quadratic alpha calculation
    bool use_linear_alpha = false;    // Use linear alpha calculation
    bool invert_transparency = false; // Invert the transparency logic

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &in_png, "Input PNG file (grayscale)")
    LONG_STRING_PARAM("out", &out_png, "Output PNG file")

    LONG_PARAM_GROUP("Color settings", NULL)
    LONG_INT_PARAM("transparent-below", &thres_transparency_below, "Threshold for transparent pixels (default: 0)")
    LONG_INT_PARAM("transparent-above", &thres_transparency_above, "Threshold for transparent pixels (default: 255)")
    LONG_STRING_PARAM("rgb", &rgb, "RGB color for the monochrome image (default: #FFFFFF)")
    LONG_PARAM("linear-alpha", &use_linear_alpha, "Use linear alpha calculation (default: constant)")
    LONG_PARAM("quadratic-alpha", &use_quadratic_alpha, "Use quadratic alpha calculation (default: constant)")
    LONG_PARAM("invert-transparency", &invert_transparency, "Invert the transparency logic (default: false)")
    LONG_INT_PARAM("alpha-constant-max", &alpha_constant_max, "Maximum value for constant alpha (default: 255)")
    LONG_DOUBLE_PARAM("alpha-quadratic-coef", &alpha_quadratic_coef, "Coefficient for quadratic alpha calculation (default: 16.0)")
    LONG_DOUBLE_PARAM("alpha-linear-coef", &alpha_linear_coef, "Coefficient for linear alpha calculation (default: 1.0)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");
    if ( in_png.empty() || out_png.empty() )
        error("--in and --out must be specified");

    const char* rgb_c = rgb.c_str();
    if ( rgb_c[0] == '#' ) ++rgb_c; // remove the # sign
    uint32_t rgb_i[3];
    sscanf(rgb_c, "%02x%02x%02x", &rgb_i[0], &rgb_i[1], &rgb_i[2]); // parse RGB hex code
    
    notice("Loading %s", in_png.c_str());
    // std::vector<uint8_t> grayscale_data; // Vector to hold grayscale data
    // uint32_t width, height;
    cimg_library::CImg<uint8_t> gray(in_png.c_str());
    if ( gray.spectrum() != 1 ) {
        error("Input PNG is not grayscale");
        return 1;
    }
    uint32_t width = gray.width();
    uint32_t height = gray.height();
    notice("Image dimensions: %dx%d", width, height);

    // if (!load_grayscale_with_fpng(in_png.c_str(), grayscale_data, width, height)) {
    //     return 1;
    // }

    // if ( grayscale_data.empty() || width == 0 || height == 0 ) {
    //     error("Error: Grayscale data size = %zu, width = %d, height = %d", grayscale_data.size(), width, height);
    //     return 1;
    // }

    // Convert grayscale to RGBA
    fpng::fpng_init(); // Initialize FPNG

    std::vector<uint8_t> rgba_data;
    rgba_data.reserve(width * height * 4);
    notice("Processing individual pixels");

    uint64_t npixels = width * height;
    uint64_t nbelow = 0, nabove = 0;
    for(int32_t i=0; i < height; ++i) {
        for(int32_t j=0; j < width; ++j) {
            uint8_t gray_value = gray(j,i); // Assuming 1 byte per pixel for grayscale
            uint32_t alpha_value = alpha_constant_max; 
            if (gray_value < thres_transparency_below) {
                rgba_data.push_back(0); // R (transparent)
                rgba_data.push_back(0); // G
                rgba_data.push_back(0); // B
                rgba_data.push_back(0); // Alpha
                ++nbelow;
            }
            else if (gray_value > thres_transparency_above) {
                rgba_data.push_back(0); // R (transparent)
                rgba_data.push_back(0); // G
                rgba_data.push_back(0); // B
                rgba_data.push_back(0); // Alpha
                ++nabove;
            }
            else {            
                rgba_data.push_back((uint8_t)floor(gray_value / 255.0 * rgb_i[0])); // R
                rgba_data.push_back((uint8_t)floor(gray_value / 255.0 * rgb_i[1])); // R
                rgba_data.push_back((uint8_t)floor(gray_value / 255.0 * rgb_i[2])); // R
                if ( use_linear_alpha ) {
                    alpha_value = (int32_t)(( invert_transparency ? 255 - gray_value : gray_value ) * alpha_linear_coef);
                }
                else if ( use_quadratic_alpha ) {
                    alpha_value = (int32_t)(sqrt((double)( invert_transparency ? 255 - gray_value : gray_value )) * alpha_quadratic_coef);
                }
                else {
                    alpha_value = alpha_constant_max; // Default constant alpha
                }
                if (alpha_value < 0) alpha_value = 0;
                if (alpha_value > 255) alpha_value = 255;
                rgba_data.push_back((uint8_t)alpha_value); 
            }
        }
    }
    notice("Processed %llu pixels, %llu below / %llu above threshold", npixels, nbelow, nabove);

    notice("Saving the png file to %s", out_png.c_str());
    if (!save_rgba_with_fpng(out_png.c_str(), rgba_data, width, height)) {
        return 1;
    }
    notice("Analysis finished");
    return 0;
}