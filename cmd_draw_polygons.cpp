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
#include <limits>
#include <fstream>
//#include <gd.h>
#define cimg_display 0  // remove the need for X11 library
//#define cimg_use_png 1 // <--- Add this line
#include "cimg/CImg.h"
#include "nlohmann/json.hpp"

using namespace cimg_library;
//using jsonnlohmann::json;

struct Color { unsigned char r, g, b; };

// Function to retrieve color based on JSON attributes
// You can customize this logic to map any attribute (ID, type, value) to a color.
Color getColorFromAttributes(const nlohmann::json& attributes) {
    // Example logic: Color based on "zone_type"
    std::string type = attributes.value("zone_type", "default");

    if (type == "water")       return {0, 100, 255};   // Blue
    if (type == "vegetation")  return {34, 139, 34};   // Forest Green
    if (type == "industrial")  return {169, 169, 169}; // Grey
    if (type == "residential") return {255, 165, 0};   // Orange
    
    return {255, 255, 255}; // Default White
}

void renderPolygon(CImg<unsigned char>& img, const nlohmann::json& item, double xmin, double ymin, double scale,
                    std::string clust_attr,
                    const std::map<std::string, uint32_t>& cluster2rgba,
                    bool cluster_as_rgb, float opacity) {    
    // print out the item for debugging
    //std::cout << item.dump(4) << std::endl; // Uncomment for
    if (!item.contains("geometry")) {
        error("Invalid Feature: missing 'geometry' field in the GeoJSON");
    }
    const auto &geometry  = item["geometry"];
    if ( !item.contains("properties") ) {
        error("Invalid Feature: missing 'properties' field in the GeoJSON");
    }
    const auto &properties = item["properties"];
    std::string clust = properties[clust_attr].get<std::string>();
    unsigned char color[] = {0, 0, 0}; // default black
    if ( cluster_as_rgb ) {
        // interpret cluster name as RRGGBB
        const char* hex = clust.c_str();
        if ( hex[0] == '#' ) hex++; // skip leading '#'
        if ( strlen(hex) != 6 ) {
            error("Invalid RGB cluster name: %s", clust.c_str());
        }
        int32_t r, g, b;
        sscanf(hex, "%2x%2x%2x", &r, &g, &b);
        color[0] = static_cast<unsigned char>(r);
        color[1] = static_cast<unsigned char>(g);
        color[2] = static_cast<unsigned char>(b);
    }
    else {
        auto it = cluster2rgba.find(clust);
        if ( it != cluster2rgba.end() ) {
            uint32_t rgba = it->second;
            color[0] = (rgba >> 16) & 0xFF; // R
            color[1] = (rgba >> 8) & 0xFF;  // G
            color[2] = rgba & 0xFF;         // B
        }
        else {
            // cluster not found, skip plotting
            return;
        }
    }

    if (!geometry.contains("type") || !geometry.contains("coordinates")) {
        error("Invalid geometry: missing 'type' or 'coordinates' field");
    }

    if (!geometry["type"].is_string()) {
        error("Invalid geometry: 'type' must be a string");
    }
    std::string type = geometry["type"].get<std::string>();
    const auto &coordinates = geometry["coordinates"][0];
    int32_t num_points = coordinates.size();
    CImg<double> points(num_points, 2);
    for (int32_t i = 0; i < num_points; i++) {
        double x = coordinates[i][0].get<double>();
        double y = coordinates[i][1].get<double>();
        // int32_t px = static_cast<int32_t>( (x - xmin) / scale );
        // int32_t py = static_cast<int32_t>( (y - ymin) / scale ); 
        points(i, 0) = (x - xmin) / scale;
        points(i, 1) = (y - ymin) / scale;
    }

    // C. Draw
    // opacity = 1.0f (fully opaque)
    img.draw_polygon(points, color, opacity);
    
    // Optional: Draw outline in black
    //unsigned char black[] = {0, 0, 0};
    //img.draw_polygon(points, black, 1.0f, ~0U); // ~0U indicates an outline pattern
}

// void renderPolygons(CImg<unsigned char>& img, const nlohmann::json& data, double xmin, double ymin, double scale,
//                     std::string clust_attr,
//                     const std::map<std::string, uint32_t>& cluster2rgba,
//                     bool cluster_as_rgb, float opacity) {    
//     for (const auto& item : data) {
//         // print out the item for debugging
//         std::cout << item.dump(4) << std::endl; // Uncomment for
//         if (!item.contains("geometry")) {
//             error("Invalid Feature: missing 'geometry' field in the GeoJSON");
//         }
//         const auto &geometry  = item["geometry"];
//         if ( !item.contains("properties") ) {
//             error("Invalid Feature: missing 'properties' field in the GeoJSON");
//         }
//         const auto &properties = item["properties"];
//         std::string clust = properties[clust_attr].get<std::string>();
//         unsigned char color[] = {0, 0, 0}; // default black
//         if ( cluster_as_rgb ) {
//             // interpret cluster name as RRGGBB
//             const char* hex = clust.c_str();
//             if ( hex[0] == '#' ) hex++; // skip leading '#'
//             if ( strlen(hex) != 6 ) {
//                 error("Invalid RGB cluster name: %s", clust.c_str());
//             }
//             int32_t r, g, b;
//             sscanf(hex, "%2x%2x%2x", &r, &g, &b);
//             color[0] = static_cast<unsigned char>(r);
//             color[1] = static_cast<unsigned char>(g);
//             color[2] = static_cast<unsigned char>(b);
//         }
//         else {
//             auto it = cluster2rgba.find(clust);
//             if ( it != cluster2rgba.end() ) {
//                 uint32_t rgba = it->second;
//                 color[0] = (rgba >> 16) & 0xFF; // R
//                 color[1] = (rgba >> 8) & 0xFF;  // G
//                 color[2] = rgba & 0xFF;         // B
//             }
//             else {
//                 // cluster not found, skip plotting
//                 continue;
//             }
//         }

//         if (!geometry.contains("type") || !geometry.contains("coordinates")) {
//             error("Invalid geometry: missing 'type' or 'coordinates' field");
//         }

//         if (!geometry["type"].is_string()) {
//             error("Invalid geometry: 'type' must be a string");
//         }
//         std::string type = geometry["type"].get<std::string>();
//         const auto &coordinates = geometry["coordinates"][0];
//         int32_t num_points = coordinates.size();
//         CImg<int> points(num_points, 2);
//         for (int32_t i = 0; i < num_points; i++) {
//             double x = coordinates[i][0].get<double>();
//             double y = coordinates[i][1].get<double>();
//             int32_t px = static_cast<int32_t>( (x - xmin) / scale );
//             int32_t py = static_cast<int32_t>( (y - ymin) / scale ); 
//             points(i, 0) = px;
//             points(i, 1) = py;
//         }

//         // C. Draw
//         // opacity = 1.0f (fully opaque)
//         img.draw_polygon(points, color, opacity);
        
//         // Optional: Draw outline in black
//         //unsigned char black[] = {0, 0, 0};
//         //img.draw_polygon(points, black, 1.0f, ~0U); // ~0U indicates an outline pattern
//     }
// }

/////////////////////////////////////////////////////////////////////////
// draw-polygons : Draw polygons, factors, and/or genes
////////////////////////////////////////////////////////////////////////
int32_t cmdDrawPolygons(int32_t argc, char **argv)
{
    std::string boundaryf;
    std::string range_tsvf;    
    double scale = 1.0; // 1 pixel = 1 coordinate unit
    int32_t verbose_freq = 100000; // report frequency of input reading
    std::string ullr;
    std::string outf;
    std::string clust_attr("topK");
    std::string cmapf;
    bool cluster_as_rgb = false;
    double opacity = 0.5;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("boundary", &boundaryf, "GeoJSON file to draw polygons from")
    LONG_STRING_PARAM("ullr", &ullr, "Specify upper-left and lower-right coordinates as 'ulx,uly,lrx,lry'")
    LONG_STRING_PARAM("range", &range_tsvf, "TSV file specifying the boundary box")
    LONG_STRING_PARAM("clust-attr", &clust_attr, "Cluster name to determine color mapping")
    LONG_STRING_PARAM("cmap", &cmapf, "RGB color mapping file")
    LONG_PARAM("cluster-as-rgb", &cluster_as_rgb, "Indicate that cluster names are RGB values as 'RRGGBB'")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_DOUBLE_PARAM("scale", &scale, "Scale factor: pixels per coordinate unit")
    LONG_STRING_PARAM("out", &outf, "Output file name")
    LONG_DOUBLE_PARAM("opacity", &opacity, "Opacity of the polygons (0.0 to 1.0)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started...");

    if (boundaryf.empty()) {
        error("Error: --boundary is required.");
        return 1;
    }
    if (outf.empty()) {
        error("Error: --out is required.");
        return 1;
    }
    if ( range_tsvf.empty() && ullr.empty() ) {
        error("Error: Either --range or --ullr is required.");
        return 1;
    }
    if ( ullr.empty() == false && range_tsvf.empty() == false ) {
        error("Error: Only one of --range or --ullr should be provided.");
        return 1;
    }
    double xmin = std::numeric_limits<double>::quiet_NaN();
    double ymin = std::numeric_limits<double>::quiet_NaN();
    double xmax = std::numeric_limits<double>::quiet_NaN();
    double ymax = std::numeric_limits<double>::quiet_NaN();
    int32_t width = 0, height = 0;
    if ( !ullr.empty() ) {
        std::vector<std::string> v;
        split(v, ",", ullr);
        if ( v.size() != 4 ) error("Invalid ullr format: %s", ullr.c_str());
        xmin = atof(v[0].c_str());
        ymin = atof(v[1].c_str());
        xmax = atof(v[2].c_str());
        ymax = atof(v[3].c_str());
    }
    else { // read range file
        tsv_reader range_tr(range_tsvf.c_str());
        while( range_tr.read_line() ) {
            const char* key = range_tr.str_field_at(0);
            double value = range_tr.double_field_at(1);
            if ( strncmp(key, "xmin", 4) == 0 ) xmin = value;
            else if ( strncmp(key, "ymin", 4) == 0 ) ymin = value;
            else if ( strncmp(key, "xmax", 4) == 0 ) xmax = value;
            else if ( strncmp(key, "ymax", 4) == 0 ) ymax = value;
        }
        if ( std::isnan(xmin) || std::isnan(ymin) || std::isnan(xmax) || std::isnan(ymax) ) {
            error("Error: Incomplete boundary box in %s", range_tsvf.c_str());
            return 1;
        }
    }
    width = (int32_t)(ceil((xmax - xmin + 1) / scale));
    height = (int32_t)(ceil((ymax - ymin + 1) / scale));
    notice("Bounding box: xmin = %lf, ymin = %lf, xmax = %lf, ymax = %lf", xmin, ymin, xmax, ymax);
    notice("Width = %d, Height = %d", width, height);

    // load color map if provided
    notice("Loading color map...");
    std::map<std::string, uint32_t> cluster2rgba;
    if ( cluster_as_rgb ) {
        notice("Interpreting cluster names as RGB values. Ignoring --cmap option even if provided.");
    }
    else {
        if ( cmapf.empty() ) {
            error("Error: --cmap is required unless --cluster-as-rgb is specified.");
            return 1;
        }
        tsv_reader cmap_tr(cmapf.c_str());        
        int32_t icol_name = -1;
        int32_t icol_r = -1;
        int32_t icol_g = -1;
        int32_t icol_b = -1;
        int32_t icol_hex = -1;
        //int32_t icol_a = -1;
        for(int32_t i = 0; cmap_tr.read_line(); i++) {
            if ( i == 0 ) {
                // identify Name, R G B
                for(int32_t j = 0; j < cmap_tr.nfields; j++) {
                    const char* h = cmap_tr.str_field_at(j);
                    if ( strncmp(h, "Name", 4) == 0 ) icol_name = j;
                    else if ( strncmp(h, "R", 1) == 0 ) icol_r = j;
                    else if ( strncmp(h, "G", 1) == 0 ) icol_g = j;
                    else if ( strncmp(h, "B", 1) == 0 ) icol_b = j;
                    //else if ( strncmp(h, "A", 1) == 0 ) icol_a = j;
                }
                if ( icol_name < 0 ) {
                    error("Error: Missing Name column in color map file %s", cmapf.c_str());
                    return 1;
                }
                if ( icol_hex < 0 && ( icol_r < 0 || icol_g < 0 || icol_b < 0 ) ) {
                    error("Error: Missing R, G, or B column in color map file %s", cmapf.c_str());
                    return 1;
                }
            }
            else {
                // process color map entries
                const char* name = cmap_tr.str_field_at(icol_name);
                if ( icol_hex >= 0 ) { // HEX color provided
                    const char* hexstr = cmap_tr.str_field_at(icol_hex);
                    if (hexstr[0] == '#' ) hexstr = hexstr + 1;
                    if ( strlen(hexstr) != 6 ) {
                        error("Error: Invalid HEX color %s for cluster %s", hexstr, name);
                        return 1;
                    }
                    int32_t r8, g8, b8;
                    sscanf(hexstr, "%2x%2x%2x", &r8, &g8, &b8);
                    uint8_t a8 = 255;
                    
                    uint32_t rgba = (a8 << 24) | (r8 << 16) | (g8 << 8) | b8;
                    cluster2rgba[std::string(name)] = rgba;
                }
                else { // RGB values provided
                    double r = cmap_tr.double_field_at(icol_r);
                    double g = cmap_tr.double_field_at(icol_g);
                    double b = cmap_tr.double_field_at(icol_b);
                    // normalize to 0-255
                    uint8_t r8 = (uint8_t)(255.0 * r);
                    uint8_t g8 = (uint8_t)(255.0 * g);
                    uint8_t b8 = (uint8_t)(255.0 * b);
                    //uint8_t a8 = icol_a >= 0 ? (uint8_t)(255.0 * cmap_tr.double_field_at(icol_a)) : 255;
                    uint8_t a8 = 255;
                    
                    uint32_t rgba = (a8 << 24) | (r8 << 16) | (g8 << 8) | b8;
                    cluster2rgba[std::string(name)] = rgba;
                }
            }
        }
    }

    // Create a blank image: Width 800, Height 600, Depth 1, Channels 3 (RGB)
    // Filled with value 20 (dark background)
    CImg<unsigned char> canvas(width, height, 1, 3, 0);

    // Load polygon data from GeoJSON
    notice("Loading polygon data from %s", boundaryf.c_str());
    nlohmann::json polygonData;
    std::ifstream inputFile(boundaryf);
    int32_t n_polygons = 0;
    if (!inputFile.is_open()) {
        error("Error: Could not open file %s", boundaryf.c_str());
        return 1;
    }
    while ( true ) {
        try {
            inputFile >> polygonData;
            renderPolygon(canvas, polygonData, xmin, ymin, scale, clust_attr, cluster2rgba, cluster_as_rgb, static_cast<float>(opacity));
            n_polygons++;
            if ( n_polygons % verbose_freq == 0 ) {
                notice("Processed %d polygons...", n_polygons);
            }
        }
        catch (nlohmann::json::parse_error& e) {
            // Check if we reached the end of the file
            if (inputFile.eof()) {
                break; // Exit the loop if end of file is reached
            } else {
                error("JSON parse error: %s", e.what());
                return 1;
            }
        }
    }
    inputFile.close();
        
    notice("Saving image to %s", outf.c_str());
    canvas.save(outf.c_str()); 

    notice("Analysis finished.");

    return 0;
}
