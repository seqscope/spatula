#ifndef __POLYGON__H__
#define __POLYGON__H__

#include <climits>
#include <cstdint>
#include <vector>
#include <fstream>
#include "nlohmann/json.hpp"
#include "qgenlib/qgen_error.h"

struct _point_t
{
    double x, y;
    _point_t(double x, double y) : x(x), y(y) {}
};
typedef struct _point_t point_t;

class Polygon
{
public:
    std::vector<point_t> vertices;
    void add_offset(double x, double y) {
        for(size_t i=0; i < vertices.size(); ++i) {
            vertices[i].x += x;
            vertices[i].y += y;
        }
    }
    bool contains_point(double x, double y) {
        size_t i, j;
        bool result = false;
        for (i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++)
        {
            if ((vertices[i].y > y) != (vertices[j].y > y) &&
                (x < (vertices[j].x - vertices[i].x) * (y - vertices[i].y) / (vertices[j].y - vertices[i].y) + vertices[i].x))
            {
                result = !result;
            }
        }
        return result;
    }

    bool contains_point(const point_t &p)
    {
        return contains_point(p.x, p.y);
    }
};

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

#endif // __POLYGON__H__