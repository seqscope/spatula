#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include <cmath>
#include <cfloat>
#include <ctime>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <thread>
#include <algorithm>
#include "nlohmann/json.hpp"

/////////////////////////////////////////////////////////////////////////
// tsv-add-cell-id : Annotate transcript-level TSV with the cell_id of the
//                   overlapping boundary polygon.
//
// The transcript TSV is streamed (never fully loaded), while the boundary
// polygons are held in memory with a uniform-grid spatial index that makes
// each point-in-polygon lookup examine only the handful of polygons whose
// bounding box overlaps the query point's grid cell. Boundaries are assumed
// disjoint, so a point normally hits 0 or 1 polygon; ties (a point inside
// more than one polygon) are broken by the closer polygon centroid.
////////////////////////////////////////////////////////////////////////

namespace {

// A single boundary polygon (exterior ring only).
struct CellPolygon {
    std::vector<double> xs;
    std::vector<double> ys;
    double xmin, ymin, xmax, ymax;   // bounding box
    double cx, cy;                    // centroid (mean of vertices)
    int32_t cell_id_idx;              // index into the global cell-id name table

    CellPolygon() : xmin(0), ymin(0), xmax(0), ymax(0), cx(0), cy(0), cell_id_idx(-1) {}

    // finalize bounding box and centroid once all vertices are loaded
    void finalize() {
        size_t n = xs.size();
        // drop a duplicated closing vertex when computing the centroid
        size_t n_centroid = (n >= 2 && xs[0] == xs[n - 1] && ys[0] == ys[n - 1]) ? n - 1 : n;
        xmin = ymin = DBL_MAX;
        xmax = ymax = -DBL_MAX;
        double sx = 0.0, sy = 0.0;
        for (size_t i = 0; i < n; ++i) {
            if (xs[i] < xmin) xmin = xs[i];
            if (xs[i] > xmax) xmax = xs[i];
            if (ys[i] < ymin) ymin = ys[i];
            if (ys[i] > ymax) ymax = ys[i];
        }
        for (size_t i = 0; i < n_centroid; ++i) { sx += xs[i]; sy += ys[i]; }
        if (n_centroid > 0) { cx = sx / n_centroid; cy = sy / n_centroid; }
    }

    // ray-casting point-in-polygon test, guarded by the bounding box
    inline bool contains(double x, double y) const {
        if (x < xmin || x > xmax || y < ymin || y > ymax) return false;
        bool in = false;
        size_t n = xs.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            if (((ys[i] > y) != (ys[j] > y)) &&
                (x < (xs[j] - xs[i]) * (y - ys[i]) / (ys[j] - ys[i]) + xs[i]))
                in = !in;
        }
        return in;
    }

    // squared distance from (px,py) to a line segment (ax,ay)-(bx,by)
    static inline double seg_dist2(double px, double py, double ax, double ay, double bx, double by) {
        double dx = bx - ax, dy = by - ay;
        double len2 = dx * dx + dy * dy;
        double t = len2 > 0.0 ? ((px - ax) * dx + (py - ay) * dy) / len2 : 0.0;
        if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
        double ex = px - (ax + t * dx), ey = py - (ay + t * dy);
        return ex * ex + ey * ey;
    }

    // squared distance from (x,y) to the nearest polygon edge (used for expansion)
    inline double dist2_to_boundary(double x, double y) const {
        double best = DBL_MAX;
        size_t n = xs.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++) {
            double d = seg_dist2(x, y, xs[j], ys[j], xs[i], ys[i]);
            if (d < best) best = d;
        }
        return best;
    }
};

// Uniform-grid spatial index over polygon bounding boxes. Because each polygon
// is registered in every grid cell its bounding box overlaps, the polygon that
// contains a query point is always registered in that point's own grid cell, so
// a query only needs to scan a single cell's candidate list. This is exact
// regardless of the cell size (cell size only trades memory for candidates).
class GridIndex {
public:
    double gxmin, gymin;
    double cell_size, inv_cell;
    int32_t ncols, nrows;
    std::vector<std::vector<int32_t> > buckets;

    void build(const std::vector<CellPolygon>& polys, double requested_cell) {
        gxmin = gymin = DBL_MAX;
        double gxmax = -DBL_MAX, gymax = -DBL_MAX;
        double sum_w = 0.0, sum_h = 0.0;
        for (const CellPolygon& p : polys) {
            if (p.xmin < gxmin) gxmin = p.xmin;
            if (p.ymin < gymin) gymin = p.ymin;
            if (p.xmax > gxmax) gxmax = p.xmax;
            if (p.ymax > gymax) gymax = p.ymax;
            sum_w += (p.xmax - p.xmin);
            sum_h += (p.ymax - p.ymin);
        }
        if (polys.empty()) { gxmin = gymin = 0; gxmax = gymax = 1; }

        // pick a cell size ~ the average polygon extent unless the user set one
        cell_size = requested_cell;
        if (cell_size <= 0.0) {
            double avg = polys.empty() ? 1.0 : (sum_w + sum_h) / (2.0 * polys.size());
            cell_size = avg > 0.0 ? avg : 1.0;
        }

        // cap the total number of grid cells to bound memory
        const int64_t MAX_CELLS = 32000000LL;
        for (;;) {
            ncols = (int32_t)((gxmax - gxmin) / cell_size) + 1;
            nrows = (int32_t)((gymax - gymin) / cell_size) + 1;
            if (ncols < 1) ncols = 1;
            if (nrows < 1) nrows = 1;
            if ((int64_t)ncols * nrows <= MAX_CELLS) break;
            cell_size *= 2.0;
        }
        inv_cell = 1.0 / cell_size;

        buckets.assign((size_t)ncols * nrows, std::vector<int32_t>());
        for (int32_t pi = 0; pi < (int32_t)polys.size(); ++pi) {
            const CellPolygon& p = polys[pi];
            int32_t c0 = clamp_col((int32_t)((p.xmin - gxmin) * inv_cell));
            int32_t c1 = clamp_col((int32_t)((p.xmax - gxmin) * inv_cell));
            int32_t r0 = clamp_row((int32_t)((p.ymin - gymin) * inv_cell));
            int32_t r1 = clamp_row((int32_t)((p.ymax - gymin) * inv_cell));
            for (int32_t r = r0; r <= r1; ++r)
                for (int32_t c = c0; c <= c1; ++c)
                    buckets[(size_t)r * ncols + c].push_back(pi);
        }
    }

    // Returns the cell-id index of the assigned polygon, or -1 if unassigned.
    // A point strictly inside a polygon is assigned by containment (ties broken
    // by the closer centroid). Otherwise, if expand > 0, the point is assigned
    // to the polygon whose boundary edge is nearest, provided it lies within
    // `expand` of it -- i.e. an exact outward buffer of the polygons by `expand`.
    inline int32_t query(double x, double y, const std::vector<CellPolygon>& polys, double expand) const {
        int32_t col = (int32_t)std::floor((x - gxmin) * inv_cell);
        int32_t row = (int32_t)std::floor((y - gymin) * inv_cell);

        // containment: only the point's own grid cell can hold a container
        if (col >= 0 && col < ncols && row >= 0 && row < nrows) {
            const std::vector<int32_t>& b = buckets[(size_t)row * ncols + col];
            int32_t best = -1;
            double best_d = DBL_MAX;
            for (size_t k = 0; k < b.size(); ++k) {
                const CellPolygon& p = polys[b[k]];
                if (p.contains(x, y)) {
                    double dx = x - p.cx, dy = y - p.cy;
                    double d = dx * dx + dy * dy;
                    if (d < best_d) { best_d = d; best = p.cell_id_idx; }
                }
            }
            if (best >= 0) return best;
        }

        if (expand <= 0.0) return -1;

        // expansion: nearest polygon edge within `expand`, searched over the
        // ceil(expand/cell) neighborhood of grid cells (exact, since a polygon
        // within `expand` has its bounding box registered in one of those cells)
        double expand2 = expand * expand;
        int32_t rad = (int32_t)std::ceil(expand * inv_cell);
        int32_t c0 = col - rad, c1 = col + rad, r0 = row - rad, r1 = row + rad;
        if (c0 < 0) c0 = 0;  if (c1 >= ncols) c1 = ncols - 1;
        if (r0 < 0) r0 = 0;  if (r1 >= nrows) r1 = nrows - 1;
        int32_t best = -1;
        double best_d = DBL_MAX;
        for (int32_t r = r0; r <= r1; ++r) {
            for (int32_t c = c0; c <= c1; ++c) {
                const std::vector<int32_t>& b = buckets[(size_t)r * ncols + c];
                for (size_t k = 0; k < b.size(); ++k) {
                    const CellPolygon& p = polys[b[k]];
                    if (x < p.xmin - expand || x > p.xmax + expand ||
                        y < p.ymin - expand || y > p.ymax + expand) continue;
                    double d = p.dist2_to_boundary(x, y);
                    if (d <= expand2 && d < best_d) { best_d = d; best = p.cell_id_idx; }
                }
            }
        }
        return best;
    }

private:
    inline int32_t clamp_col(int32_t c) const { return c < 0 ? 0 : (c >= ncols ? ncols - 1 : c); }
    inline int32_t clamp_row(int32_t r) const { return r < 0 ? 0 : (r >= nrows ? nrows - 1 : r); }
};

// Maps cell-id strings to a compact index, so each polygon stores an int and
// the streaming stage can print names without repeated string hashing.
class CellIdTable {
public:
    std::vector<std::string> names;
    std::unordered_map<std::string, int32_t> name2idx;

    int32_t get_or_add(const std::string& name) {
        std::unordered_map<std::string, int32_t>::iterator it = name2idx.find(name);
        if (it != name2idx.end()) return it->second;
        int32_t idx = (int32_t)names.size();
        names.push_back(name);
        name2idx[name] = idx;
        return idx;
    }
};

// Convert a GeoJSON scalar (string or number) property into a string cell id.
std::string json_scalar_to_string(const nlohmann::json& v) {
    if (v.is_string()) return v.get<std::string>();
    if (v.is_number_integer()) return std::to_string(v.get<int64_t>());
    if (v.is_number_unsigned()) return std::to_string(v.get<uint64_t>());
    if (v.is_number_float()) {
        char buf[64];
        snprintf(buf, sizeof(buf), "%g", v.get<double>());
        return std::string(buf);
    }
    return v.dump();
}

// Append one polygon ring (array of [x,y]) as a CellPolygon.
void add_ring(const nlohmann::json& ring, int32_t cell_id_idx, std::vector<CellPolygon>& polys) {
    if (!ring.is_array() || ring.size() < 3) return;
    CellPolygon poly;
    poly.cell_id_idx = cell_id_idx;
    poly.xs.reserve(ring.size());
    poly.ys.reserve(ring.size());
    for (size_t i = 0; i < ring.size(); ++i) {
        const nlohmann::json& pt = ring[i];
        if (!pt.is_array() || pt.size() < 2 || !pt[0].is_number() || !pt[1].is_number())
            error("Invalid GeoJSON position: expected [x, y] numbers");
        poly.xs.push_back(pt[0].get<double>());
        poly.ys.push_back(pt[1].get<double>());
    }
    poly.finalize();
    polys.push_back(std::move(poly));
}

// Add all polygons contained in a single GeoJSON Feature.
void add_feature(const nlohmann::json& feature, const std::string& cell_id_attr,
                 CellIdTable& cell_ids, std::vector<CellPolygon>& polys) {
    if (!feature.contains("geometry")) error("Invalid Feature: missing 'geometry'");
    const nlohmann::json& geometry = feature["geometry"];
    if (!geometry.contains("type") || !geometry.contains("coordinates"))
        error("Invalid geometry: missing 'type' or 'coordinates'");
    if (!feature.contains("properties") || !feature["properties"].contains(cell_id_attr))
        error("Invalid Feature: missing property '%s'", cell_id_attr.c_str());

    std::string type = geometry["type"].get<std::string>();
    int32_t cid = cell_ids.get_or_add(json_scalar_to_string(feature["properties"][cell_id_attr]));
    const nlohmann::json& coords = geometry["coordinates"];

    if (type == "Polygon") {
        if (coords.is_array() && !coords.empty())
            add_ring(coords[0], cid, polys);   // exterior ring only; holes ignored
    } else if (type == "MultiPolygon") {
        for (size_t i = 0; i < coords.size(); ++i)
            if (coords[i].is_array() && !coords[i].empty())
                add_ring(coords[i][0], cid, polys);
    } else {
        error("Unsupported geometry type '%s' (expected Polygon or MultiPolygon)", type.c_str());
    }
}

// Load polygons from a GeoJSON file. Handles both a single FeatureCollection
// object and a stream of concatenated Feature objects (one per line/record).
void load_geojson(const std::string& fn, const std::string& cell_id_attr,
                  CellIdTable& cell_ids, std::vector<CellPolygon>& polys) {
    std::ifstream in(fn.c_str());
    if (!in.is_open()) error("Cannot open GeoJSON file %s", fn.c_str());
    uint64_t n_features = 0;
    while (true) {
        nlohmann::json j;
        try {
            in >> j;
        } catch (nlohmann::json::parse_error& e) {
            if (in.eof()) break;   // clean end of stream
            error("GeoJSON parse error in %s: %s", fn.c_str(), e.what());
        }
        if (j.is_object() && j.contains("type") && j["type"] == "FeatureCollection") {
            for (const nlohmann::json& feat : j["features"]) {
                add_feature(feat, cell_id_attr, cell_ids, polys);
                if (++n_features % 100000 == 0) notice("Loaded %llu features...", n_features);
            }
        } else if (j.is_object() && j.contains("type") && j["type"] == "Feature") {
            add_feature(j, cell_id_attr, cell_ids, polys);
            if (++n_features % 100000 == 0) notice("Loaded %llu features...", n_features);
        } else {
            error("Invalid GeoJSON %s: expected 'Feature' or 'FeatureCollection'", fn.c_str());
        }
    }
    notice("Loaded %llu features, %zu polygons from %s", n_features, polys.size(), fn.c_str());
}

// Strip a single pair of surrounding double quotes, if present.
inline const char* unquote(const char* s, std::string& buf) {
    size_t len = strlen(s);
    if (len >= 2 && s[0] == '"' && s[len - 1] == '"') {
        buf.assign(s + 1, len - 2);
        return buf.c_str();
    }
    return s;
}

// Load polygons from a CSV(.gz) with cell_id / vertex_x / vertex_y columns.
// Vertices belonging to the same cell_id (in file order) form one polygon.
void load_csv(const std::string& fn, const std::string& col_cell_id,
              const std::string& col_x, const std::string& col_y,
              CellIdTable& cell_ids, std::vector<CellPolygon>& polys) {
    tsv_reader tr(fn.c_str());
    tr.delimiter = ',';
    if (tr.read_line() <= 0) error("Empty CSV boundary file %s", fn.c_str());

    int32_t i_cid = -1, i_x = -1, i_y = -1;
    std::string qb;
    for (int32_t i = 0; i < tr.nfields; ++i) {
        const char* h = unquote(tr.str_field_at(i), qb);
        if (col_cell_id == h) i_cid = i;
        else if (col_x == h) i_x = i;
        else if (col_y == h) i_y = i;
    }
    if (i_cid < 0 || i_x < 0 || i_y < 0)
        error("CSV %s must contain columns '%s', '%s', '%s'", fn.c_str(),
              col_cell_id.c_str(), col_x.c_str(), col_y.c_str());

    // map cell-id string -> index of the polygon currently being built
    std::unordered_map<std::string, int32_t> cid2poly;
    std::string cbuf, xbuf, ybuf;
    uint64_t n_rows = 0;
    while (tr.read_line() > 0) {
        std::string cid = unquote(tr.str_field_at(i_cid), cbuf);
        double x = atof(unquote(tr.str_field_at(i_x), xbuf));
        double y = atof(unquote(tr.str_field_at(i_y), ybuf));
        std::unordered_map<std::string, int32_t>::iterator it = cid2poly.find(cid);
        int32_t pi;
        if (it == cid2poly.end()) {
            pi = (int32_t)polys.size();
            cid2poly[cid] = pi;
            polys.push_back(CellPolygon());
            polys.back().cell_id_idx = cell_ids.get_or_add(cid);
        } else {
            pi = it->second;
        }
        polys[pi].xs.push_back(x);
        polys[pi].ys.push_back(y);
        if (++n_rows % 1000000 == 0) notice("Loaded %llu vertices...", n_rows);
    }
    tr.close();
    for (size_t i = 0; i < polys.size(); ++i) polys[i].finalize();
    notice("Loaded %llu vertices, %zu polygons from %s", n_rows, polys.size(), fn.c_str());
}

// Locate the start of the xi-th and yi-th tab-delimited fields in one scan.
// strtod stops at the trailing tab, so the raw line needs no modification.
inline bool extract_xy(const std::string& line, int32_t xi, int32_t yi, double& x, double& y) {
    const char* s = line.c_str();
    const char* px = NULL;
    const char* py = NULL;
    int32_t f = 0;
    const char* start = s;
    for (const char* p = s;; ++p) {
        if (*p == '\t' || *p == '\0') {
            if (f == xi) px = start;
            if (f == yi) py = start;
            if (*p == '\0') break;
            ++f;
            start = p + 1;
        }
    }
    if (px == NULL || py == NULL) return false;
    x = strtod(px, NULL);
    y = strtod(py, NULL);
    return true;
}

} // namespace

int32_t cmdTsvAddCellId(int32_t argc, char **argv) {
    std::string tsvf;
    std::string geojsonf;
    std::string csvf;
    std::string outf;
    std::string colname_x("X");
    std::string colname_y("Y");
    std::string out_colname_cell_id("cell_id");
    std::string unassigned_cell_id("UNASSIGNED");
    std::string cell_id_attr("cell_id");
    std::string csv_col_cell_id("cell_id");
    std::string csv_col_x("vertex_x");
    std::string csv_col_y("vertex_y");
    double bucket_um = 0.0;   // 0 = auto (average polygon extent)
    double expand_um = 0.0;   // outward buffer distance for boundaries (0 = off)
    int32_t nthreads = (int32_t)std::thread::hardware_concurrency();
    int32_t batch_size = 500000;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Transcript TSV Options", NULL)
    LONG_STRING_PARAM("tsv", &tsvf, "Transcript-level TSV(.gz) file with X/Y coordinates (streamed)")
    LONG_STRING_PARAM("colname-x", &colname_x, "Column name for X coordinate in the TSV")
    LONG_STRING_PARAM("colname-y", &colname_y, "Column name for Y coordinate in the TSV")

    LONG_PARAM_GROUP("Input Boundary Options (one required)", NULL)
    LONG_STRING_PARAM("boundaries-geojson", &geojsonf, "GeoJSON file with cell boundary polygons")
    LONG_STRING_PARAM("boundaries-csv", &csvf, "CSV(.gz) with cell_id, vertex_x, vertex_y columns")
    LONG_STRING_PARAM("cell-id-attr", &cell_id_attr, "GeoJSON property name holding the cell ID")
    LONG_STRING_PARAM("csv-colname-cell-id", &csv_col_cell_id, "CSV column name for cell ID")
    LONG_STRING_PARAM("csv-colname-x", &csv_col_x, "CSV column name for vertex X")
    LONG_STRING_PARAM("csv-colname-y", &csv_col_y, "CSV column name for vertex Y")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output TSV(.gz) file with the added cell-id column")
    LONG_STRING_PARAM("out-colname-cell-id", &out_colname_cell_id, "Column name for the added cell ID")
    LONG_STRING_PARAM("unassigned-cell-id", &unassigned_cell_id, "Value for transcripts inside no polygon")

    LONG_PARAM_GROUP("Boundary Expansion Options", NULL)
    LONG_DOUBLE_PARAM("expand-um", &expand_um, "Expand cell boundaries outward by this distance; a point outside all polygons is assigned to the nearest polygon within this distance (0 = off)")

    LONG_PARAM_GROUP("Performance Options", NULL)
    LONG_INT_PARAM("threads", &nthreads, "Number of threads for point-in-polygon assignment")
    LONG_INT_PARAM("batch-size", &batch_size, "Number of TSV lines processed per batch")
    LONG_DOUBLE_PARAM("bucket-um", &bucket_um, "Grid cell size for the spatial index (0 = auto)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (tsvf.empty()) error("--tsv is required");
    if (outf.empty()) error("--out is required");
    if (geojsonf.empty() == csvf.empty())
        error("Exactly one of --boundaries-geojson or --boundaries-csv must be specified");
    if (nthreads < 1) nthreads = 1;
    if (batch_size < 1) batch_size = 1;
    if (expand_um < 0.0) error("--expand-um must be non-negative");
    if (expand_um > 0.0)
        notice("Boundaries will be expanded outward by %.4g (nearest-edge assignment)", expand_um);

    // ---- load boundary polygons into memory -------------------------------
    CellIdTable cell_ids;
    std::vector<CellPolygon> polys;
    if (!geojsonf.empty()) {
        notice("Loading polygons from GeoJSON %s", geojsonf.c_str());
        load_geojson(geojsonf, cell_id_attr, cell_ids, polys);
    } else {
        notice("Loading polygons from CSV %s", csvf.c_str());
        load_csv(csvf, csv_col_cell_id, csv_col_x, csv_col_y, cell_ids, polys);
    }
    if (polys.empty()) error("No polygons were loaded from the boundary file");

    // ---- build the spatial index ------------------------------------------
    notice("Building spatial index over %zu polygons...", polys.size());
    GridIndex grid;
    grid.build(polys, bucket_um);
    notice("Spatial index: %d x %d grid, cell size %.4g", grid.ncols, grid.nrows, grid.cell_size);

    // ---- open input and output; parse header ------------------------------
    htsFile* hp = hts_open(tsvf.c_str(), "r");
    if (hp == NULL) error("Cannot open TSV file %s", tsvf.c_str());
    htsFile* wf = hts_open(outf.c_str(),
                           outf.size() >= 3 && outf.compare(outf.size() - 3, 3, ".gz") == 0 ? "wz" : "w");
    if (wf == NULL) error("Cannot open output file %s for writing", outf.c_str());

    kstring_t str = {0, 0, NULL};
    if (hts_getline(hp, KS_SEP_LINE, &str) <= 0) error("TSV file %s is empty", tsvf.c_str());

    int32_t x_idx = -1, y_idx = -1;
    {   // tokenize the header line by tab
        int32_t f = 0;
        char* start = str.s;
        for (char* p = str.s;; ++p) {
            if (*p == '\t' || *p == '\0') {
                char saved = *p;
                *p = '\0';
                if (colname_x == start) x_idx = f;
                if (colname_y == start) y_idx = f;
                if (out_colname_cell_id == start)
                    error("Column '%s' already exists in %s; use --out-colname-cell-id",
                          out_colname_cell_id.c_str(), tsvf.c_str());
                *p = saved;
                if (saved == '\0') break;
                ++f;
                start = p + 1;
            }
        }
    }
    if (x_idx < 0 || y_idx < 0)
        error("Columns '%s' and/or '%s' not found in the TSV header", colname_x.c_str(), colname_y.c_str());

    hprintf(wf, "%s\t%s\n", str.s, out_colname_cell_id.c_str());

    // ---- stream the TSV in batches, assign cell ids in parallel -----------
    std::vector<std::string> lines;
    lines.reserve(batch_size);
    std::vector<int32_t> results;
    uint64_t n_lines = 0, n_assigned = 0;

    while (true) {
        // read a batch of raw lines
        lines.clear();
        for (int32_t i = 0; i < batch_size; ++i) {
            if (hts_getline(hp, KS_SEP_LINE, &str) <= 0) break;
            lines.push_back(std::string(str.s, str.l));
        }
        if (lines.empty()) break;

        int32_t n = (int32_t)lines.size();
        results.assign(n, -1);

        // partition the batch across worker threads
        int32_t nt = nthreads < n ? nthreads : n;
        auto worker = [&](int32_t t) {
            int32_t lo = (int64_t)n * t / nt;
            int32_t hi = (int64_t)n * (t + 1) / nt;
            for (int32_t i = lo; i < hi; ++i) {
                double x, y;
                if (extract_xy(lines[i], x_idx, y_idx, x, y))
                    results[i] = grid.query(x, y, polys, expand_um);
            }
        };
        if (nt == 1) {
            worker(0);
        } else {
            std::vector<std::thread> pool;
            pool.reserve(nt);
            for (int32_t t = 0; t < nt; ++t) pool.emplace_back(worker, t);
            for (std::thread& th : pool) th.join();
        }

        // write results in input order
        for (int32_t i = 0; i < n; ++i) {
            const std::string& name = results[i] < 0 ? unassigned_cell_id : cell_ids.names[results[i]];
            if (results[i] >= 0) ++n_assigned;
            hprintf(wf, "%s\t%s\n", lines[i].c_str(), name.c_str());
        }
        n_lines += n;
        notice("Processed %llu transcripts (%llu assigned to a cell)...", n_lines, n_assigned);
    }

    hts_close(wf);
    hts_close(hp);
    if (str.s) free(str.s);

    notice("Finished: %llu transcripts, %llu assigned to a cell (%.2f%%)",
           n_lines, n_assigned, n_lines > 0 ? 100.0 * n_assigned / n_lines : 0.0);
    notice("Analysis finished");
    return 0;
}
