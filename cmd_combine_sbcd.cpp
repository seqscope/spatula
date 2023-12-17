#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge.h"
#include "seq_utils.h"
#include "file_utils.h"
#include <ctime>
#include <cmath>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

struct _tile_info_t
{
    int32_t row;
    int32_t col;
    uint64_t xmin;
    uint64_t xmax;
    uint64_t ymin;
    uint64_t ymax;
    uint64_t x_offset;
    uint64_t y_offset;
    std::string sbcdf;

    _tile_info_t() : row(0), col(0), xmin(0), xmax(0), ymin(0), ymax(0), x_offset(0), y_offset(0) {}
};
typedef struct _tile_info_t tile_info_t;

///////////////////////////////////////////////////////////////////////////////////////////////
// combine-sge : Combine multiple SDGE tiles with global coordinates based on specified layouts
///////////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdCombineSBCD(int32_t argc, char **argv)
{
    std::string layoutf;             // Layout file, each containing [lane] [tile] and [row]/[col] as columns
    std::string offsetf;             // Offset file, each containing [lane] [tile] and [x_offset]/[y_offset] as columns
    std::string sbcddir;             // Directory containing spatial barcode files
    std::string manifestf;           // Manifest file containing the list of spatial barcode files
    std::string outdir;              // output directory containing merged sbcd with global coordinates in nm scale
    int32_t match_len = 27;          // length of HDMI spatial barcodes to be considered for matching
    double pixel_to_nm = 34.78;      // pixel to nm conversion factor (37.5 for Seq-Scope Hi-Seq, 34.78 for Seq-Scope NovaSeq)
    int32_t max_dup_allowed = 5;     // maximum number of duplicates allowed for each spatial barcode. If this is 1, duplicates are not allowed
    double max_dup_dist_nm = 10000.; // maximum distance allowed for duplicates in nm scale
    double rowgap = 0.0;             // additional gap between rows (proportional to the height of a tile)
    double colgap = 0.0;             // additional gap between columns (proportional to the width of a tile)
    bool write_all = false;          // write all spatial barcodes to the output file

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("layout", &layoutf, "Layout file, each containing [lane] [tile] and [row]/[col] as columns")
    LONG_STRING_PARAM("offset", &offsetf, "Offset file, each containing [lane] [tile] and [row]/[col] as columns")
    LONG_STRING_PARAM("sbcd", &sbcddir, "Directory containing spatial barcode files")
    LONG_STRING_PARAM("manifest", &manifestf, "Manifest file containing the list of spatial barcode files")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outdir, "Output spatial barcode file after merging")
    LONG_PARAM("write-all", &write_all, "Write all spatial barcodes to the output file, including duplicated and filtered reads")

    LONG_PARAM_GROUP("Options for coordinate conversion", NULL)
    LONG_DOUBLE_PARAM("pixel-to-nm", &pixel_to_nm, "Pixel to nm conversion factor (37.5 for Seq-Scope)")
    LONG_DOUBLE_PARAM("rowgap", &rowgap, "Additional gap between rows (proportional to the height of a tile)")
    LONG_DOUBLE_PARAM("colgap", &colgap, "Additional gap between columns (proportional to the width of a tile)")

    LONG_PARAM_GROUP("Options for duplicate filtering", NULL)
    LONG_INT_PARAM("match-len", &match_len, "Length of HDMI spatial barcode to be considered for a match")
    LONG_INT_PARAM("max-dup", &max_dup_allowed, "Maximum number of duplicates allowed for each spatial barcode. If this is 1, duplicates are not allowed")
    LONG_DOUBLE_PARAM("max-dup-dist-nm", &max_dup_dist_nm, "Maximum distance allowed for duplicates in nm scale")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (sbcddir.empty() || manifestf.empty() || outdir.empty())
    {
        error("Missing required options --sbcd, --manifest, --out");
    }
    if ( layoutf.empty() == offsetf.empty() ) {
        error("Only one option should be set between --layout and --offset");
    }

    dataframe_t manifest_df(manifestf.c_str()); // open the manifest file

    // layout file is expected to have [lane] [tile] [row] [col], and optionally [rowshift] [colshift]
    // in an advanced layout, layout file can have [lane] [tile] [x_offset] [y_offset]
    // manifest file is expected to have [id] [filepath] [barcodes] [matches] [mismatches] [xmin] [xmax] [ymin] [ymax]

    std::map<std::string, tile_info_t*> tile_info_map;

    // scan the manifest file and store the max/min coordinates
    uint64_t max_xdiff = 0, max_ydiff = 0;
    uint64_t min_xdiff = UINT64_MAX, min_ydiff = UINT64_MAX;
    char buf_id[255], buf_filepath[65535];
    int32_t i_xmin = manifest_df.get_colidx("xmin");
    int32_t i_xmax = manifest_df.get_colidx("xmax");
    int32_t i_ymin = manifest_df.get_colidx("ymin");
    int32_t i_ymax = manifest_df.get_colidx("ymax");
    if ( i_xmin < 0 || i_xmax < 0 || i_ymin < 0 || i_ymax < 0 ) {
        error("Manifest file %s does not have xmin, xmax, ymin, ymax columns", manifestf.c_str());
    }
    for (int32_t i = 0; i < manifest_df.nrows; ++i)
    {
        std::string sid(manifest_df.get_str_elem(i, 0));
        snprintf(buf_filepath, 65535, "%s/%s", sbcddir.c_str(), manifest_df.get_str_elem(i, 1).c_str());
        uint64_t xmin = manifest_df.get_uint64_elem(i, i_xmin);
        uint64_t xmax = manifest_df.get_uint64_elem(i, i_xmax);
        uint64_t ymin = manifest_df.get_uint64_elem(i, i_ymin);
        uint64_t ymax = manifest_df.get_uint64_elem(i, i_ymax);
        if (xmax - xmin + 1 > max_xdiff)
            max_xdiff = xmax - xmin + 1;
        if (ymax - ymin + 1 > max_ydiff)
            max_ydiff = ymax - ymin + 1;
        if (xmax - xmin + 1 < min_xdiff)
            min_xdiff = xmax - xmin + 1;
        if (ymax - ymin + 1 < min_ydiff)
            min_ydiff = ymax - ymin + 1;

        // set attribute for the tile
        if ( tile_info_map.find(sid) == tile_info_map.end() ) {
            tile_info_t* pti = new tile_info_t;
            pti->xmin = xmin;
            pti->xmax = xmax;
            pti->ymin = ymin;
            pti->ymax = ymax;
            pti->sbcdf = buf_filepath;
            tile_info_map[sid] = pti;
        }
        else {
            error("Duplicate lane_tile %s exists in the manifest file %s", buf_id, manifestf.c_str());
        }
    }

    // add a column to indicate the full path
    int32_t icol = manifest_df.add_empty_column("fullpath");
    int32_t jcol = manifest_df.get_colidx("filepath");
    for(int32_t i=0; i < manifest_df.nrows; ++i) {
        manifest_df.set_str_elem((sbcddir + "/" + manifest_df.get_str_elem(i, jcol)).c_str(), i, icol);
    }


    notice("Finished reading the manifest file across %d tiles", manifest_df.nrows);
    notice("max_xdiff = %llu, min_xdiff = %llu, difference = %llu", max_xdiff, min_xdiff, max_xdiff-min_xdiff);
    notice("max_ydiff = %llu, min_ydiff = %llu, difference = %llu", max_ydiff, min_ydiff, max_ydiff-min_ydiff);

    std::map<std::string, tile_info_t*>::iterator it;
    if ( !layoutf.empty() ) {
        // First, read the layout file and find out the rules for coordinate conversion
        dataframe_t layout_df(layoutf.c_str());     // open the layout file

        int32_t i_id = layout_df.get_colidx("id");
        int32_t i_lane = layout_df.get_colidx("lane");
        int32_t i_tile = layout_df.get_colidx("tile");
        int32_t i_row = layout_df.get_colidx("row");
        int32_t i_col = layout_df.get_colidx("col");
        int32_t i_rowshift = layout_df.get_colidx("rowshift");
        int32_t i_colshift = layout_df.get_colidx("colshift");

        if ( i_row < 0 || i_col < 0 ) 
            error("[row] and [col] are required but missing in the layout file %s", layoutf.c_str());

        // find out the maximum row value
        int32_t max_row = 0;
        for (int32_t i = 0; i < layout_df.nrows; ++i) {
            int32_t row = layout_df.get_int_elem(i, i_row);
            max_row = row > max_row ? row : max_row;
        }

        // Read the layout file to determine the offsets for each tile
        // calculate offset in the following way
        // x_offset = height * ((1 + rowgap) * (row - 1) + rowshift)
        // y_offset = width *  ((1 + colgap) * (col - 1) + colshift)
        for (int32_t i = 0; i < layout_df.nrows; ++i)
        {
            tile_info_t* pti = NULL;
            //if ( layout_df.has_column("id") ) {
            if ( i_id >= 0 ) {
                it = tile_info_map.find(layout_df.get_str_elem(i, i_id));
                if ( it == tile_info_map.end() ) {
                    error("Tile %s does not exist in the layout file %s", layout_df.get_str_elem(i, i_id).c_str(), layoutf.c_str());
                }
                pti = it->second;
            }
            else {
                //if ( !( layout_df.has_column("lane") && layout_df.has_column("tile") ) ) {
                if ( i_lane < 0 || i_tile < 0 ) {
                    error("[id] or [lane]/[tile] column is required in the layout file %s", layoutf.c_str());
                }
                snprintf(buf_id, 255, "%d_%d", layout_df.get_int_elem(i, i_lane), layout_df.get_int_elem(i, i_tile));
                it = tile_info_map.find(buf_id);
                if ( it == tile_info_map.end() ) {
                    error("Tile %s does not exist in the layout file %s", buf_id, layoutf.c_str());
                }
                pti = it->second;
            }

            int32_t row = layout_df.get_int_elem(i, i_row);
            int32_t col = layout_df.get_int_elem(i, i_col);
            double rowshift = i_rowshift < 0 ? 0. : layout_df.get_double_elem(i, i_rowshift);
            double colshift = i_colshift < 0 ? 0. : layout_df.get_double_elem(i, i_colshift);
            //uint64_t x_offset = (uint64_t)(max_xdiff * ( (rowgap + 1.0) * (max_row - row) + rowshift ));
            uint64_t x_offset = (uint64_t)(max_xdiff * ( (rowgap + 1.0) * (row - 1) + rowshift ));
            uint64_t y_offset = (uint64_t)(max_ydiff * ( (colgap + 1.0) * (col - 1) + colshift ));
            pti->x_offset = x_offset; // now we can just add 
            pti->y_offset = y_offset;
        }
    }
    else {
        // First, read the offset file and find out the rules for coordinate conversion
        dataframe_t offset_df(offsetf.c_str());     // open the offset file

        int32_t i_id = offset_df.get_colidx("id");
        int32_t i_lane = offset_df.get_colidx("lane");
        int32_t i_tile = offset_df.get_colidx("tile");
        int32_t i_xoffset = offset_df.get_colidx("x_offset");
        int32_t i_yoffset = offset_df.get_colidx("y_offset");

        if ( ( i_xoffset < 0 ) || ( i_yoffset < 0 ) )
            error("[x_offset] and [y_offset] are required but missing in the offset file %s", offsetf.c_str());

        for (int32_t i = 0; i < offset_df.nrows; ++i)
        {
            tile_info_t* pti = NULL;
            if ( i_id >= 0 ) {
                it = tile_info_map.find(offset_df.get_str_elem(i, i_id));
                if ( it == tile_info_map.end() ) {
                    error("Tile %s does not exist in the manifest file %s", offset_df.get_str_elem(i, i_id).c_str(), manifestf.c_str());
                }
                pti = it->second;
            }
            else {
                if ( i_lane < 0 || i_tile < 0 ) {
                    error("[id] or [lane]/[tile] column is required in the offset file %s", offsetf.c_str());
                }
                snprintf(buf_id, 255, "%d_%d", offset_df.get_int_elem(i, i_lane), offset_df.get_int_elem(i, i_tile));
                it = tile_info_map.find(buf_id);
                if ( it == tile_info_map.end() ) {
                    error("Tile %s does not exist in the manifest file %s", buf_id, manifestf.c_str());
                }
                pti = it->second;
            }

            pti->x_offset = (uint64_t)offset_df.get_uint64_elem(i, i_xoffset);
            pti->y_offset = (uint64_t)offset_df.get_uint64_elem(i, i_yoffset);
        }        
    }   

    notice("Finished calculating offsets from input files");

    // create the output directory first
    if (makePath(outdir))
    {
        notice("Successfully created the directory %s", outdir.c_str());
    }
    else
    {
        notice("The directory %s already exists", outdir.c_str());
    }

    // write a new sbcd file
    htsFile* wh_sbcd = hts_open( (outdir + "/1_1.sbcds.sorted.tsv.gz").c_str(), "wz" );
    htsFile* wh_manifest = hts_open( (outdir + "/manifest.tsv").c_str(), "w" );
    htsFile* wh_dupstat = hts_open( (outdir + "/dupstats.tsv.gz").c_str(), "wz" );
    htsFile* wh_dups = write_all ? hts_open( (outdir + "/1_1.dups.sorted.tsv.gz").c_str(), "wz" ) : NULL;
    htsFile* wh_filt = write_all ? hts_open( (outdir + "/1_1.filtered.sorted.tsv.gz").c_str(), "wz" ) : NULL;

    hprintf(wh_manifest, "id\tfilepath\tbarcodes\ttmatches\tmismatches\txmin\txmax\tymin\tymax\n");
    hprintf(wh_dupstat,  "ndups\tmax_dist_nm\n");

    // open all tiles
    std::vector<std::string> tiles;
    std::vector<tsv_reader*> bcdfs;
    open_tiles(manifest_df, tiles, bcdfs);
    int32_t ntiles = manifest_df.nrows;
    std::vector<uint64_t> tseqs(ntiles);   // sequences at each tile

    // read first lines from each tile
    uint64_t min_tseq = UINT64_MAX;
    std::vector<int32_t> imins_tseq;
    for (int32_t i = 0; i < ntiles; ++i)
    {
        tseqs[i] = seq2nt5(bcdfs[i]->str_field_at(0), match_len);
        if ( min_tseq > tseqs[i] ) {
            min_tseq = tseqs[i];
            imins_tseq.clear();
            imins_tseq.push_back(i);
        }
        else if ( min_tseq == tseqs[i] ) {
            imins_tseq.push_back(i);
        }
    }

    std::vector<sbcd_rec_t> sbcd_recs;
    uint64_t npass = 0, nfilt = 0, ufilt = 0, ndups = 0, udups = 0;
    uint64_t min_gx = UINT64_MAX, max_gx = 0, min_gy = UINT64_MAX, max_gy = 0;
    while( min_tseq < UINT64_MAX ) {
        // collect all ties;
        for(int32_t i=0; i < (int32_t)imins_tseq.size(); ++i) {
            // insert the current element to the vector
            int32_t j = imins_tseq[i];
            while ( tseqs[j] == min_tseq ) {
                sbcd_recs.emplace_back(tseqs[j],    // nid
                    bcdfs[j]->str_field_at(0),      // strid
                    bcdfs[j]->int_field_at(1),      // lane
                    bcdfs[j]->int_field_at(2),      // tile
                    bcdfs[j]->uint64_field_at(3),   // px
                    bcdfs[j]->uint64_field_at(4),   // py
                    bcdfs[j]->int_field_at(5));     // mismatch
                // read the next line
                if (bcdfs[j]->read_line() == 0) // EOF reached
                    tseqs[j] = UINT64_MAX;
                else    // valid barcode
                    tseqs[j] = seq2nt5(bcdfs[j]->str_field_at(0), match_len);
            }
        }

        // perform duplicate filtering
        if ( sbcd_recs.size() > max_dup_allowed ) {
            // filter the duplicate reads
            if ( write_all ) {
                for(int32_t i=0; i < (int32_t)sbcd_recs.size(); ++i) {
                    sbcd_recs[i].hprint_sbcd(wh_filt);
                }
            }
            hprintf(wh_dupstat, "%zu\tNA\n", sbcd_recs.size());
            nfilt += sbcd_recs.size();
            ++ufilt;
        }
        else {
            // convert the current coordinates into global coordinates
            for(int32_t i=0; i < sbcd_recs.size(); ++i) {
                sbcd_rec_t& rec = sbcd_recs[i];
                snprintf(buf_id, 255, "%llu_%llu", rec.lane, rec.tile);
                it = tile_info_map.find(buf_id);
                if ( it == tile_info_map.end() ) {
                    error("Tile ID %s does not exist in the manifest file %s", buf_id, manifestf.c_str());
                }
                tile_info_t* pti = tile_info_map[buf_id];
                //uint64_t gx = (uint64_t)((rec.px + pti->x_offset - pti->xmin) * pixel_to_nm);
                uint64_t gx = (uint64_t)((pti->x_offset + pti->xmax - rec.px) * pixel_to_nm);
                uint64_t gy = (uint64_t)((rec.py + pti->y_offset - pti->ymin) * pixel_to_nm);
                if ( gx < min_gx ) min_gx = gx;
                if ( gx > max_gx ) max_gx = gx;
                if ( gy < min_gy ) min_gy = gy;
                if ( gy > max_gy ) max_gy = gy;
                rec.px = gy; // swap x/y to make it compatible with visualization
                rec.py = gx; 
                rec.lane = 1;
                rec.tile = 1;
            }

            // calculate the maximum distance between the barcodes
            double max_dist_nm = 0;
            for(int32_t i=0; i < sbcd_recs.size(); ++i) {
                for(int32_t j=0; j < i; ++j) {
                    double xdiff = (double)sbcd_recs[i].px - (double)sbcd_recs[j].px;
                    double ydiff = (double)sbcd_recs[i].py - (double)sbcd_recs[j].py;
                    double dist_nm = sqrt( xdiff*xdiff + ydiff*ydiff );
                    if ( dist_nm > max_dist_nm ) max_dist_nm = dist_nm;
                }
            }

            if ( sbcd_recs.size() > 1 )
                hprintf(wh_dupstat, "%zu\t%llu\n", sbcd_recs.size(), (uint64_t)max_dist_nm);

            // decide whether to keep the reads or not
            if ( max_dist_nm > max_dup_dist_nm ) {
                // duplicate reads are too far away
                if ( write_all ) {
                    for(int32_t i=0; i < (int32_t)sbcd_recs.size(); ++i) {
                        sbcd_recs[i].lane = (uint64_t) sbcd_recs.size();
                        sbcd_recs[i].tile = (uint64_t) max_dist_nm;
                        sbcd_recs[i].hprint_sbcd(wh_filt);
                    }
                }
                nfilt += sbcd_recs.size();
                ++ufilt;
            }
            else {  // the reads are close enough to be considered as optical duplicates
                // in this case, write an arbitrary (first) read to the output
                sbcd_recs[0].hprint_sbcd(wh_sbcd);
                ++npass;
                if ( sbcd_recs.size() > 1 ) {
                    ++udups;
                    ndups += (sbcd_recs.size()-1);
                    if ( write_all ) {
                        // write the rest of reads to the duplicate output
                        for(int32_t i=1; i < sbcd_recs.size(); ++i) {
                            sbcd_recs[i].lane = (uint64_t) sbcd_recs.size();
                            sbcd_recs[i].tile = (uint64_t) max_dist_nm;
                            sbcd_recs[i].hprint_sbcd(wh_dups);
                        }
                    }
                }
            }
        }
        // determine the next min_tseq
        imins_tseq.clear();
        min_tseq = UINT64_MAX;
        for (int32_t i = 0; i < ntiles; ++i)
        {
            if ( min_tseq > tseqs[i] ) {
                min_tseq = tseqs[i];
                imins_tseq.clear();
                imins_tseq.push_back(i);
            }
            else if ( min_tseq == tseqs[i] ) {
                imins_tseq.push_back(i);
            }
        }
        sbcd_recs.clear();

        if ( npass % 1000000 == 0 ) {
            notice("Processed %llu records; %llu passed, %llu (%llu unique) filtered, %llu (%llu unique) duplicates", npass+nfilt+ndups, npass, nfilt, ufilt, ndups, udups);
        }
    }
    notice("Finished processing %llu records; %llu passed, %llu (%llu unique) filtered, %llu (%llu unique) duplicates", npass+nfilt+ndups, npass, nfilt, ufilt, ndups, udups);

    hts_close(wh_sbcd);
    hts_close(wh_dupstat);
    if ( write_all ) {
        hts_close(wh_filt);
        hts_close(wh_dups);
    }

    hprintf(wh_manifest, "1_1\t1_1.sbcds.sorted.tsv.gz\t%llu\t%llu\t0\t%llu\t%llu\t%llu\t%llu\n", npass, npass, min_gy, max_gy, min_gx, max_gx); // note x/y are swapped
    hts_close(wh_manifest);

    notice("Analysis finished");

    return 0;
}
