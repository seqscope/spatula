#include "polygon.h"

int32_t load_polygons_from_geojson(const char *filename, std::vector<Polygon> &polygons)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }

    nlohmann::json json;
    file >> json;

    //notice("foo");
    auto &feature = json["features"][0];

    //notice("bar");

    Polygon polygon;
    auto &poly_coords = feature["geometry"]["coordinates"];

    //notice("baz");

    for (auto &coordinates : poly_coords)
    {
        Polygon polygon;
        for (auto &coordinate : coordinates[0])
        {
            double x = coordinate[0];
            double y = coordinate[1];
            polygon.vertices.push_back(point_t(x, y));
        }
        polygons.push_back(polygon);
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
