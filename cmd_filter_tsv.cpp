#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "file_utils.h"
#include "polygon.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
// filter-tsv : Filter TSV based on spatial coordinates or gene list
////////////////////////////////////////////////////////////////////////////////
int32_t cmdFilterTSV(int32_t argc, char **argv)
{
    std::string in_tsv;
    std::string out_prefix;

    std::string tsv_suffix("transcripts.tsv.gz");
    std::string ftr_suffix("features.tsv.gz");
    std::string minmax_suffix("minmax.tsv");

    // output column names for tsv files
    std::string colname_X("X");
    std::string colname_Y("Y");
    std::string colname_gene("gene");
    std::string colname_gid("gene_id");
    std::string colname_cnt("gn");

    // filtering option
    double in_xmin = -std::numeric_limits<double>::infinity();
    double in_xmax = std::numeric_limits<double>::infinity();
    double in_ymin = -std::numeric_limits<double>::infinity();
    double in_ymax = std::numeric_limits<double>::infinity();

    std::string filt_genef;
    std::string filt_gidf;
    std::string geojsonf;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("in-tsv", &in_tsv, "Input (unsorted) TSV file")
    LONG_STRING_PARAM("out-prefix", &out_prefix, "Prefix of output files")
    LONG_STRING_PARAM("tsv-suffix", &tsv_suffix, "Suffix for output TSV file")
    LONG_STRING_PARAM("ftr-suffix", &ftr_suffix, "Suffix for output feature file")
    LONG_STRING_PARAM("minmax-suffix", &ftr_suffix, "Suffix for output minmax file")

    LONG_PARAM_GROUP("Column Names", NULL)
    LONG_STRING_PARAM("colname-x", &colname_X, "Column name for X-axis")
    LONG_STRING_PARAM("colname-y", &colname_Y, "Column name for Y-axis")
    LONG_STRING_PARAM("colname-gene", &colname_gene, "Column name for gene name")
    LONG_STRING_PARAM("colname-gid", &colname_gid, "Column name for gene ID")

    LONG_PARAM_GROUP("Spatial Filtering options", NULL)
    LONG_DOUBLE_PARAM("xmin", &in_xmin, "Minimum x-axis value")
    LONG_DOUBLE_PARAM("xmax", &in_xmax, "Maximum x-axis value")
    LONG_DOUBLE_PARAM("ymin", &in_ymin, "Minimum y-axis value")
    LONG_DOUBLE_PARAM("ymax", &in_ymax, "Maximum y-axis value")
    LONG_STRING_PARAM("polygon", &geojsonf, "GeoJSON file for polygon-based filtering")

    LONG_PARAM_GROUP("Gene Filtering options", NULL)
    LONG_STRING_PARAM("filt-gene", &filt_genef, "Only Include gene names present in the list")
    LONG_STRING_PARAM("filt-gid", &filt_gidf, "Only Include gene ids present in the list")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( in_tsv.empty() )
        error("Input TSV, must be specified");

    notice("Analysis started");

    // read the input feature
    std::set<std::string> set_gname;
    std::set<std::string> set_gid;

    if ( !filt_genef.empty() ) {
        tsv_reader tr_gene(filt_genef.c_str());
        while( tr_gene.read_line() ) {
            set_gname.insert(tr_gene.str_field_at(0));
        }
    }

    if ( !filt_gidf.empty() ) {
        tsv_reader tr_gid(filt_gidf.c_str());
        while( tr_gid.read_line() ) {
            set_gid.insert(tr_gid.str_field_at(0));
        }
    }

    // read polygon filter
    int32_t npolygons = 0;
    std::vector<Polygon> polygons;
    if (!geojsonf.empty())
    {
        npolygons = load_polygons_from_geojson(geojsonf.c_str(), polygons);
    }

    // process input TSV file and separate into batches based on the Y-axis
    tsv_reader tr(in_tsv.c_str());
    if ( !tr.read_line() )
        error("Failed to read the header line from %s", in_tsv.c_str());

    htsFile* wh = hts_open((out_prefix + tsv_suffix).c_str(), "wz");
    if ( wh == NULL )
        error("Failed to open %s%s for writing", out_prefix.c_str(), tsv_suffix.c_str());

    // identify the header column
    int32_t icol_x = -1;
    int32_t icol_y = -1;
    int32_t icol_gene = -1;
    int32_t icol_gid = -1;
    int32_t icol_cnt = -1;
    std::vector<std::string> v_colnames;
    for(int32_t i=0; i < tr.nfields; ++i) {
        v_colnames.push_back(tr.str_field_at(i));
        if ( colname_X.compare(v_colnames[i]) == 0 ) { // found the column
            if ( icol_x < 0 ) {
                icol_x = i;
            }
            else {
                error("Column %s is not unique in the header line", colname_X.c_str());
            }
        }
        if ( colname_Y.compare(v_colnames[i]) == 0 ) { // found the column
            if ( icol_y < 0 ) {
                icol_y = i;
            }
            else {
                error("Column %s is not unique in the header line", colname_Y.c_str());
            }
        }
        if ( colname_gene.compare(v_colnames[i]) == 0 ) { // found the column
            if ( icol_gene < 0 ) {
                icol_gene = i;
            }
            else {
                error("Column %s is not unique in the header line", colname_gene.c_str());
            }
        }
        if ( colname_gid.compare(v_colnames[i]) == 0 ) { // found the column
            if ( icol_gid < 0 ) {
                icol_gid = i;
            }
            else {
                error("Column %s is not unique in the header line", colname_gid.c_str());
            }
        }
        if ( colname_cnt.compare(v_colnames[i]) == 0 ) { // found the column
            if ( icol_cnt < 0 ) {
                icol_cnt = i;
            }
            else {
                error("Column %s is not unique in the header line", colname_cnt.c_str());
            }
        }

        for(int32_t i=0; i < tr.nfields; ++i) {
            hprintf(wh, "%s", tr.str_field_at(i));
            if ( i < tr.nfields - 1 )
                hprintf(wh, "\t");
        }
        hprintf(wh, "\n");
    }

    if ( icol_x < 0 )
        error("Cannot find the name of X column %s", colname_X.c_str());
    if ( icol_y < 0 )
        error("Cannot find the name of Y column %s", colname_Y.c_str());
    if ( ( icol_gene < 0 ) && ( !filt_genef.empty() ) )
        error("Cannot find the name of gene column %s", colname_gene.c_str());
    if ( ( icol_gid < 0 ) && ( !filt_gidf.empty() ) )
        error("Cannot find the name of gene ID column %s", colname_gid.c_str());
    if ( ( icol_gid < 0 ) && ( icol_gene < 0 ) ) 
        error("Either gene (%s) or gene ID (%s) column must be present in TSV", colname_gene.c_str(), colname_gid.c_str());
    if ( icol_cnt < 0 )
        error("Cannot find the name of count column %s", colname_cnt.c_str());

    // write the input TSV file into multiple small batches
    std::map<std::string, uint64_t> gene2cnt;
    double out_xmin = std::numeric_limits<double>::infinity();
    double out_xmax = -std::numeric_limits<double>::infinity();
    double out_ymin = std::numeric_limits<double>::infinity();
    double out_ymax = -std::numeric_limits<double>::infinity();
    uint64_t nlines = 0;

    bool has_spatial_filter = 
        ( in_xmin > -std::numeric_limits<double>::infinity() ) || 
        ( in_xmax <  std::numeric_limits<double>::infinity() ) || 
        ( in_ymin > -std::numeric_limits<double>::infinity() ) || 
        ( in_ymax < std::numeric_limits<double>::infinity() ) ||
        ( npolygons > 0 );
    
    bool has_gene_filter = !filt_genef.empty() || !filt_gidf.empty();

    while( tr.read_line() ) {
        bool pass = true;
        double x = tr.double_field_at(icol_x);
        double y = tr.double_field_at(icol_y);

        if ( has_spatial_filter ) {
            if ( ( x < in_xmin ) || ( x > in_xmax ) || ( y < in_ymin ) || ( y > in_ymax ) ) {
                pass = false;
            }
            else if ( npolygons > 0 ) {
                pass = false;
                for (int32_t i = 0; i < npolygons; ++i)
                {
                    if (polygons[i].contains_point(x, y))
                    {
                        pass = true;
                        break;
                    }
                }
            }
        }
        if ( has_gene_filter ) {
            if ( !filt_genef.empty() ) {
                if ( set_gname.find(tr.str_field_at(icol_gene)) == set_gname.end() ) {
                    pass = false;
                }
            }
            if ( !filt_gidf.empty() ) {
                if ( set_gid.find(tr.str_field_at(icol_gid)) == set_gid.end() ) {
                    pass = false;
                }
            }
        }

        if ( pass ) {
            // update gene count
            std::string gkey;
            if ( ( icol_gene >= 0 ) && ( icol_gid >= 0 ) ) {
                gkey = std::string(tr.str_field_at(icol_gene)) + ":" + tr.str_field_at(icol_gid);
            }
            else if ( icol_gene >= 0 ) {
                gkey = std::string(tr.str_field_at(icol_gene)) + ":" + tr.str_field_at(icol_gene);
            }
            else if ( icol_gid >= 0 ) {
                gkey = std::string(tr.str_field_at(icol_gid)) + ":" + tr.str_field_at(icol_gid);
            }
            gene2cnt[gkey] += tr.int_field_at(icol_cnt);

            // update the bounding box
            out_xmin = out_xmin < x ? out_xmin : x;
            out_xmax = out_xmax > x ? out_xmax : x;
            out_ymin = out_ymin < y ? out_ymin : y;
            out_ymax = out_ymax > y ? out_ymax : y;

            // write the line
            for(int32_t i=0; i < tr.nfields; ++i) {
                hprintf(wh, "%s", tr.str_field_at(i));
                if ( i < tr.nfields - 1 )
                    hprintf(wh, "\t");
            }
            hprintf(wh, "\n");
        }
    } 
    hts_close(wh);

    notice("Analysis finished");

    return 0;
}
