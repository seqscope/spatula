#ifndef __POLYGON__H__
#define __POLYGON__H__

#include <climits>
#include <cstdint>
#include <vector>
#include <fstream>
//#include <iostream>
#include "nlohmann/json.hpp"
#include "qgenlib/qgen_error.h"

struct _point_t
{
    double x, y;
    _point_t(double x, double y) : x(x), y(y) {}
};
typedef struct _point_t point_t;

class Rectangle
{
public:
    point_t p_min, p_max;

    Rectangle(double xmin, double ymin, double xmax, double ymax) : p_min(xmin, ymin), p_max(xmax, ymax) {}

    inline void add_point(double x, double y) {
        if (x < p_min.x) p_min.x = x;
        if (x > p_max.x) p_max.x = x;
        if (y < p_min.y) p_min.y = y;
        if (y > p_max.y) p_max.y = y;
    }

    inline bool contains_point(double x, double y) const
    {
        return (x >= p_min.x && x <= p_max.x && y >= p_min.y && y <= p_max.y);
    }

    inline bool contains_point(const point_t &p) const
    {
        return contains_point(p.x, p.y);
    }

    inline bool contains_rectangle(const Rectangle &r) const
    {
        return (contains_point(r.p_min) && contains_point(r.p_max));
    }

    inline bool intersects_rectangle(const Rectangle &r) const
    {
        return (contains_point(r.p_min) || contains_point(r.p_max) || r.contains_point(p_min) || r.contains_point(p_max));
    }
};

class Polygon
{
public:
    std::vector<point_t> vertices;
    inline void add_offset(double x, double y)
    {
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            vertices[i].x += x;
            vertices[i].y += y;
        }
    }
    inline bool contains_point(double x, double y)
    {
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

    inline bool contains_point(const point_t &p)
    {
        return contains_point(p.x, p.y);
    }

    Rectangle get_bounding_box()
    {
        double xmin = std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::min();
        double ymax = std::numeric_limits<double>::min();
        for (auto &vertex : vertices)
        {
            //notice("Vertex: %lf, %lf, (%lg %lg %lg %lg)", vertex.x, vertex.y, xmin, ymin, xmax, ymax);
            if (vertex.x < xmin)
                xmin = vertex.x;
            if (vertex.x > xmax)
                xmax = vertex.x;
            if (vertex.y < ymin)
                ymin = vertex.y;
            if (vertex.y > ymax)
                ymax = vertex.y;
        }
        return Rectangle(xmin, ymin, xmax, ymax);
    }
};

int32_t add_feature_to_polygons(const nlohmann::json &feature, std::vector<Polygon> &polygons);
int32_t load_polygons_from_geojson(const char *filename, std::vector<Polygon> &polygons);
bool polygons_contain_point(std::vector<Polygon> &polygons, double x, double y);

#endif // __POLYGON__H__