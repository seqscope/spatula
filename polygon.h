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

int32_t load_polygons_from_geojson(const char *filename, std::vector<Polygon> &polygons);
bool polygons_contain_point(std::vector<Polygon> &polygons, double x, double y);

#endif // __POLYGON__H__