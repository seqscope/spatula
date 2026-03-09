#include "polygon.h"

int32_t add_feature_to_polygons(const nlohmann::json &feature, std::vector<Polygon> &polygons) {
    if (!feature.contains("geometry")) {
        error("Invalid Feature: missing 'geometry' field in the GeoJSON");
    }

    const auto &geometry = feature["geometry"];
    if (!geometry.contains("type") || !geometry.contains("coordinates")) {
        error("Invalid geometry: missing 'type' or 'coordinates' field");
    }

    if (!geometry["type"].is_string()) {
        error("Invalid geometry: 'type' must be a string");
    }
    std::string type = geometry["type"].get<std::string>();

    const auto &coordinates = geometry["coordinates"];

    auto add_single_polygon = [&](const nlohmann::json &coords) {
        if (!coords.is_array() || coords.empty()) {
            error("Invalid Polygon: 'coordinates' must be a non-empty array of linear rings");
        }

        // We only support the exterior ring (index 0); holes are not supported.
        uint64_t n_vertices = 0;
        for (size_t ringIndex = 0; ringIndex < coords.size(); ++ringIndex) {
            if (ringIndex > 0) {
                error("Holes are not supported");
            }
            const auto &ring = coords[ringIndex];
            if (!ring.is_array() || ring.empty()) {
                error("Invalid Polygon ring: must be a non-empty array of positions");
            }

            Polygon polygon;
            polygon.vertices.reserve(ring.size());
            for (size_t pointIndex = 0; pointIndex < ring.size(); ++pointIndex) {
                const auto &pt = ring[pointIndex];
                if (!pt.is_array() || pt.size() < 2 || !pt[0].is_number() || !pt[1].is_number()) {
                    error("Invalid position: expected [x, y] numbers");
                }
                polygon.vertices.push_back(point_t(pt[0].get<double>(), pt[1].get<double>()));
                n_vertices++;
            }
            polygons.push_back(std::move(polygon));
        }
        notice("Total of %llu vertices added for a Polygon", n_vertices);
    };

    if (type == "Polygon") {
        add_single_polygon(coordinates);
    } else if (type == "MultiPolygon") {
        if (!coordinates.is_array() || coordinates.empty()) {
            error("Invalid MultiPolygon: 'coordinates' must be a non-empty array of Polygons");
        }
        for (const auto &polyCoords : coordinates) {
            add_single_polygon(polyCoords);
        }
    } else {
        error("Geometry type is not Polygon or MultiPolygon");
    }

    return static_cast<int32_t>(polygons.size());
}


// inline int32_t add_feature_to_polygons(const nlohmann::json &feature, std::vector<Polygon> &polygons)
// {
//     if (!feature.contains("geometry")) {
//         error("Invalid Feature: missing 'geometry' field in the GeoJSON");
//     }

//     auto &geometry = feature["geometry"];
    
//     if (!geometry.contains("type") || !geometry.contains("coordinates")) {
//         error("Invalid geometry: missing 'type' or 'coordinates' field");
//     }

//     if (geometry["type"] != "Polygon") {
//         error("Geometry type is not Polygon");
//     }

//     const auto& coordinates = geometry["coordinates"];

//     // Iterate through each ring (exterior and holes)
//     Polygon polygon;
//     for (size_t ringIndex = 0; ringIndex < coordinates.size(); ++ringIndex) {
//         if (ringIndex > 0) {
//             error("Holes are not supported");
//         }
            
//         // Iterate through points in the ring
//         for (size_t pointIndex = 0; pointIndex < coordinates[ringIndex].size(); ++pointIndex) {
//             const auto& point = coordinates[ringIndex][pointIndex];
//             polygon.vertices.push_back(point_t(point[0], point[1]));
//         }
//     }
//     polygons.push_back(polygon);
//     return (int32_t)polygons.size();
// }

int32_t load_polygons_from_geojson(const char *filename, std::vector<Polygon> &polygons)
{
    try {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file");
        }

        nlohmann::json json;
        file >> json;

        // search for "type" or "geometry" in the JSON file
        if ( !json.contains("type") && json.contains("geometry") ) { // assume that the JSON is not a geojson but contains geojson as "geometry"
            json = json["geometry"];
        }        

        // now the type should exist
        if (!json.contains("type")) {
             error("Invalid GeoJSON %s: missing 'type' field", filename);
        }

        if (json["type"] != "FeatureCollection" && json["type"] != "Feature") {
            error("Invalid GeoJSON %s: expected 'FeatureCollection' or 'Feature' type", filename);
        }

        uint64_t n_features = 0;
        if (json["type"] == "FeatureCollection") {
            for (const auto& feature : json["features"]) {
                add_feature_to_polygons(feature, polygons);
                n_features++;
            }
        } else if (json["type"] == "Feature") {
            add_feature_to_polygons(json, polygons);
            n_features++;
        } else {
            error("Invalid GeoJSON %s: Cannot find 'Feature' of 'FeatureCollection' type", filename);
        }
        notice("Loaded %llu features, total of %zu polygons", n_features, polygons.size());
    }
    catch (const std::exception &e) {
        error("Error loading GeoJSON %s: %s", filename, e.what());
    }
    return (int32_t)polygons.size();
}


bool polygons_contain_point(std::vector<Polygon> &polygons, double x, double y)
{
    for (auto &polygon : polygons)
    {
        if (polygon.contains_point(x, y))
        {
            return true;
        }
    }
    return false;
}
