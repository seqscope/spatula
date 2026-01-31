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
#include "nlohmann/json.hpp"

void updateBoundingBox(const nlohmann::json& item, double& xmin, double& ymin, double& xmax, double& ymax) {
    if (!item.contains("geometry")) {
        error("Invalid Feature: missing 'geometry' field in the GeoJSON");
    }
    const auto &geometry  = item["geometry"];
    if (!geometry.contains("type") || !geometry.contains("coordinates")) {
        error("Invalid geometry: missing 'type' or 'coordinates' field");
    }
    if (!geometry["type"].is_string()) {
        error("Invalid geometry: 'type' must be a string");
    }
    std::string type = geometry["type"].get<std::string>();
    const auto &coordinates = geometry["coordinates"][0];
    for (const auto& coord : coordinates) {
        double x = coord[0].get<double>();
        double y = coord[1].get<double>();
        if (x < xmin) xmin = x;
        if (y < ymin) ymin = y;
        if (x > xmax) xmax = x;
        if (y > ymax) ymax = y;
    }
}


/////////////////////////////////////////////////////////////////////////
// draw-polygons : Draw polygons, factors, and/or genes
////////////////////////////////////////////////////////////////////////
int32_t cmdPixels2Cells(int32_t argc, char **argv)
{
    std::string boundaryf;
    std::string pixelf;
    std::string range_tsvf;
    std::string ullr;
    std::string outf;
    std::string cell_id_attr("cell_id");
    std::string colname_pixel_x("X");
    std::string colname_pixel_y("Y");
    std::string out_colname_cell_id("cell_id");
    std::string unassigned_cell_id("UNASSIGNED");
    double um_per_grid = 1.0; // grid size in microns

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input GeoJSON options", NULL)
    LONG_STRING_PARAM("boundary", &boundaryf, "GeoJSON file containing cell IDs and polygons")
    LONG_STRING_PARAM("cell-id-attr", &cell_id_attr, "Attribute name for cell IDs in the GeoJSON properties")

    LONG_PARAM_GROUP("Input Pixel TSV options", NULL)
    LONG_STRING_PARAM("pixel", &pixelf, "Pixel-level TSV data containing individual molecule positions")
    LONG_STRING_PARAM("colname-x", &colname_pixel_x, "Column name for pixel X coordinate")
    LONG_STRING_PARAM("colname-y", &colname_pixel_y, "Column name for pixel Y coordinate")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output TSV file name with cell assignments")
    LONG_STRING_PARAM("out-colname-cell-id", &out_colname_cell_id, "Column name for cell ID in the output TSV")
    LONG_STRING_PARAM("unassigned-cell-id", &unassigned_cell_id, "Cell ID to assign for unassigned pixels")

    LONG_PARAM_GROUP("Other Input Options", NULL)
    LONG_DOUBLE_PARAM("um-per-grid", &um_per_grid, "Grid size in microns (default: 1.0)")
    LONG_STRING_PARAM("range", &range_tsvf, "TSV file specifying xmin, ymin, xmax, ymax for the area of interest. If both --range and --ullr is absent, the boundary file will be read twice to determine the bounding box.")
    LONG_STRING_PARAM("ullr", &ullr, "Specify upper-left and lower-right coordinates as 'ulx,uly,lrx,lry'. If both --range and --ullr is absent, the boundary file will be read twice to determine the bounding box.")
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

    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymax = std::numeric_limits<double>::lowest();
    int32_t width = 0, height = 0;
    if ( !ullr.empty() ) {
        notice("Using --ullr %s to set bounding box.", ullr.c_str());
        std::vector<std::string> v;
        split(v, ",", ullr);
        if ( v.size() != 4 ) error("Invalid ullr format: %s", ullr.c_str());
        xmin = atof(v[0].c_str());
        ymin = atof(v[1].c_str());
        xmax = atof(v[2].c_str());
        ymax = atof(v[3].c_str());
    }
    else if ( !range_tsvf.empty() ) { // read range file
        notice("Using --range %s to set bounding box.", range_tsvf.c_str());
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
    else {
        notice("Determining bounding box from boundary file %s", boundaryf.c_str());
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
                updateBoundingBox(polygonData, xmin, ymin, xmax, ymax);
                n_polygons++;
                if ( n_polygons % 100000 == 0 ) {
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
        notice("Finished processed %d polygons", n_polygons);
        //notice("Determined bounding box: xmin=%lf, ymin=%lf, xmax=%lf, ymax=%lf", xmin, ymin, xmax, ymax);
        inputFile.close();
    }
    width = (int32_t)(ceil((xmax - xmin + 1) / um_per_grid));
    height = (int32_t)(ceil((ymax - ymin + 1) / um_per_grid));
    notice("Bounding box: xmin = %lf, ymin = %lf, xmax = %lf, ymax = %lf", xmin, ymin, xmax, ymax);
    notice("Width = %d, Height = %d", width, height);

    // create a canvas and assign cell IDs to pixels
    std::vector<std::string> cell_ids;
    cell_ids.push_back(unassigned_cell_id); // index 0 for unassigned
    uint32_t** grids = (uint32_t**)malloc(sizeof(uint32_t*) * width);
    for(int32_t i = 0; i < width; i++) {
        grids[i] = (uint32_t*)calloc(height, sizeof(uint32_t)); // set everything to 0
    }

    // Load polygon data from GeoJSON to fill in grids with cell IDs
    notice("Loading polygon data from %s", boundaryf.c_str());
    nlohmann::json polygonData;
    std::ifstream inputFile(boundaryf);
    int32_t n_polygons = 0;
    uint64_t n_inconsistent_grids = 0;
    uint64_t n_assigned_grids = 0;
    if (!inputFile.is_open()) {
        error("Error: Could not open file %s", boundaryf.c_str());
        return 1;
    }
    while ( true ) {
        try {
            inputFile >> polygonData;
            if (!polygonData.contains("geometry")) {
                error("Invalid Feature: missing 'geometry' field in the GeoJSON");
            }
            const auto &geometry  = polygonData["geometry"];
            if (!geometry.contains("type") || !geometry.contains("coordinates")) {
                error("Invalid geometry: missing 'type' or 'coordinates' field");
            }
            if (!geometry["type"].is_string()) {
                error("Invalid geometry: 'type' must be a string");
            }
            std::string type = geometry["type"].get<std::string>();
            const auto &coordinates = geometry["coordinates"][0];
            const auto &properties = polygonData["properties"];
            std::string cell_id = properties[cell_id_attr].get<std::string>();
            uint32_t int_cell_id = static_cast<uint32_t>(cell_ids.size());
            cell_ids.push_back(cell_id);
            // fill grids with cell IDs
            int32_t num_points = coordinates.size();

            // fill in grids with cell IDs
            // the current bounding box is xmin, ymin, xmax, ymax
            // the grid resolution is um_per_grid
            // uint32_t** grids are allocated already
            if ( num_points >= 3) {
                double min_y =  std::numeric_limits<double>::infinity();
                double max_y = -std::numeric_limits<double>::infinity();
                std::vector<double> xs;
                std::vector<double> ys;
                for(int32_t i = 0; i < num_points; i++) {
                    const auto& p = coordinates[i];
                    xs.push_back( (p[0].get<double>() - xmin ) / um_per_grid );
                    ys.push_back( (p[1].get<double>() - ymin ) / um_per_grid );
                    min_y = std::min(min_y, ys.back());
                    max_y = std::max(max_y, ys.back());
                }
                // ensure closed ring
                if (xs.front() != xs.back() || ys.front() != ys.back()) {
                    xs.push_back(xs.front());
                    ys.push_back(ys.front());
                    ++num_points;
                }
                
                // scan-line rasterization
                int32_t start_y = std::max(0, static_cast<int>(std::floor(min_y)));
                int32_t end_y = std::min(height-1, static_cast<int>(std::floor(max_y)));
                for (int y = start_y; y <= end_y; ++y) {
                    std::vector<double> nodeX;
                    size_t j = xs.size() - 1; // Previous vertex index

                    // Find intersections with the current scanline 'y'
                    for (size_t i = 0; i < xs.size(); ++i) {
                        double Yi = ys[i];
                        double Yj = ys[j];
                        double Xi = xs[i];
                        double Xj = xs[j];

                        // Check if edge crosses the scanline
                        if ((Yi < y && Yj >= y) || (Yj < y && Yi >= y)) {
                            // Calculate X intersection
                            double intersect_x = Xi + (y - Yi) / (Yj - Yi) * (Xj - Xi);
                            nodeX.push_back(intersect_x);
                        }
                        j = i;
                    }

                    // Sort intersections from left to right
                    std::sort(nodeX.begin(), nodeX.end());

                    // Fill pixels between pairs of intersections
                    for (size_t i = 0; i < nodeX.size(); i += 2) {
                        if (i + 1 >= nodeX.size()) break;

                        // Determine integer pixel span
                        int start_x = static_cast<int>(std::ceil(nodeX[i]));
                        int end_x   = static_cast<int>(std::floor(nodeX[i + 1]));

                        // Clamp X to grid width
                        if (start_x < 0) start_x = 0;
                        if (end_x >= width) end_x = width - 1;

                        // Fill the span with the ID
                        for (int x = start_x; x <= end_x; ++x) {
                            // assign cell ID to the grid
                            if ( grids[x][y] == 0 ) { // no ID assigned
                                grids[x][y] = int_cell_id;
                                n_assigned_grids++;
                            }
                            else if ( grids[x][y] == int_cell_id ) {
                                // already assigned, do nothing
                            }
                            else if ( grids[x][y] == UINT32_MAX ) {
                                // inconsistently assigned cluster, leave it as is
                            }
                            else {
                                // inconsistently assigned cluster, mark it as MAX_UINT32
                                grids[x][y] = UINT32_MAX;
                                n_inconsistent_grids++;
                            }
                        }
                    }
                }
            }
            n_polygons++;
            if ( n_polygons % 100000 == 0 ) {
                notice("Processed %d polygons... Assigned %llu grids (%.2lf%% of total area) and identified %llu inconsistent grids", n_polygons, n_assigned_grids, static_cast<double>(n_assigned_grids) / (width * height) * 100, n_inconsistent_grids);
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
    notice("Finished processing %d polygons... Assigned %llu grids (%.2lf%% of total area) and identified %llu inconsistent grids", n_polygons, n_assigned_grids, static_cast<double>(n_assigned_grids) / (width * height) * 100, n_inconsistent_grids);
    inputFile.close();    

    // read the input TSV file and assign cell IDs to the grid
    tsv_reader tr(pixelf.c_str());
    htsFile* wf = hts_open(outf.c_str(), outf.compare(outf.size() - 3, 3, ".gz", 3) == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Failed to open output file %s", outf.c_str());
        return 1;
    }
    bool is_header = true;
    int32_t x_idx = -1;
    int32_t y_idx = -1;
    uint64_t n_lines = 0;
    while ( tr.read_line() ) {
        if ( is_header ) {
            is_header = false;
            // write header
            for (int32_t i = 0; i < tr.nfields; ++i) {
                const char* s = tr.str_field_at(i);
                if ( colname_pixel_x.compare(s) == 0 ) {
                    x_idx = i;
                }
                else if ( colname_pixel_y.compare(s) == 0 ) {
                    y_idx = i;
                }
                if ( out_colname_cell_id.compare(s) == 0 ) {
                    error("Column %s already exists in the input TSV file. Please use --out-colname-cell-id to specify a different output column name.", out_colname_cell_id.c_str());
                    return 1;
                }
                hprintf(wf, "%s\t", s);
            }
            hprintf(wf, "%s\n", out_colname_cell_id.c_str());
            if ( x_idx == -1 || y_idx == -1 ) {
                error("Column %s or %s not found in the input TSV file", colname_pixel_x.c_str(), colname_pixel_y.c_str());
                return 1;
            }
        }
        else {
            // write the output line
            for (int32_t i = 0; i < tr.nfields; ++i) {
                const char* s = tr.str_field_at(i);
                hprintf(wf, "%s\t", s);
            }
            double x = tr.double_field_at(x_idx);
            double y = tr.double_field_at(y_idx);
            int32_t grid_x = static_cast<int32_t>(std::floor((x - xmin) / um_per_grid));
            int32_t grid_y = static_cast<int32_t>(std::floor((y - ymin) / um_per_grid));
            uint32_t int_cell_id = 0;
            if ( grid_x >= 0 && grid_x < width && grid_y >= 0 && grid_y < height ) {
                int_cell_id = grids[grid_x][grid_y];
                if ( int_cell_id == UINT32_MAX ) {
                    int_cell_id = 0;
                }
            }
            if ( int_cell_id >= cell_ids.size() ) {
                error("Cell ID %u is out of bounds. This should not happen.", int_cell_id);
                return 1;
            }
            hprintf(wf, "%s\n", cell_ids[int_cell_id].c_str());
        }
        ++n_lines;
        if ( n_lines % 1000000 == 0 ) {
            notice("Processed %llu lines...", n_lines);
        }
    }
    notice("Finished processing %llu lines.", n_lines);
    hts_close(wf);
    tr.close();

    // free memory
    for (int i = 0; i < width; i++) {
        free(grids[i]);
    }
    free(grids);

    notice("Analysis finished.");

    return 0;
}
