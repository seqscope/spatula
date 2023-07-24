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


/////////////////////////////////////////////////////////////////////////
// eval-dups-sbcds : Evaluate duplicates in spatial barcodes
////////////////////////////////////////////////////////////////////////
int32_t cmdEvalDupsSBCD(int32_t argc, char **argv)
{
    std::string bcddir;
    std::string outprefix;
    int32_t match_len = 27;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("sbcd", &bcddir, "Spatial barcode dictionary generated from 'build-sbcds' command")
    LONG_INT_PARAM("match-len", &match_len, "Length of HDMI spatial barcodes to require perfect matches")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outprefix, "Output prefix (index.tsv, matches.tsv.gz)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (bcddir.empty() || outprefix.empty())
    {
        error("Missing required options --sbcd, --out");
    }

    // make spatial barcode directory end with '/'
    if (bcddir[bcddir.size() - 1] != '/')
        bcddir += "/";

    // read the manifest file
    dataframe_t df((bcddir + "manifest.tsv").c_str());
    if (df.nrows == 0)
        error("Empty dataframe %smanifest.tsv", bcddir.c_str());

    notice("Successfully read the manifest file, containing %d rows and %d columns", df.nrows, df.ncols);

    // add a column to indicate the full path
    int32_t icol = df.add_empty_column("fullpath");
    int32_t jcol = df.get_colidx("filepath");
    for (int32_t i = 0; i < df.nrows; ++i)
    {
        df.set_str_elem((bcddir + df.get_str_elem(i, jcol)).c_str(), i, icol);
    }

    std::vector<std::string> tiles;
    std::vector<tsv_reader *> bcdfs;

    open_tiles(df, tiles, bcdfs);

    int32_t ntiles = (int32_t)tiles.size();

    int32_t len = strlen(bcdfs[0]->str_field_at(0));
    if (len < match_len)
        error("HDMI length %d does not match to the parameters %d", len, match_len);

    // read the first entry in each tile, and store the minimum values and locations
    std::vector<uint64_t> tseqs(ntiles);
    std::vector<int32_t> imins;
    uint64_t nt5min = UINT64_MAX; 
    std::string seqmin;
    std::vector<uint64_t> ntotal_tiles(ntiles, 0);
    for (int32_t i = 0; i < ntiles; ++i)
    {
        tseqs[i] = read_bcdf(bcdfs[i], match_len, ntotal_tiles[i]);
        if ( nt5min > tseqs[i] ) {
            nt5min = tseqs[i];
            imins.clear();
            imins.push_back(i);
        } else if ( nt5min == tseqs[i] ) {
            imins.push_back(i);
        }
    }
    seqmin.assign(bcdfs[imins[0]]->str_field_at(0)); // assign the minimum sequence

    // output file storing individual duplicate reads
    char buf[65536];
    sprintf(buf, "%s.dups.sbcds.tsv.gz", outprefix.c_str());
    htsFile *wdups = hts_open(buf, "wz");
    hprintf(wdups, "#barcode\tntiles\tdupcount\n");

    uint64_t nuniq = 0, ndups = 0, ndups_within = 0, ndups_between = 0, ndups_uniq = 0, dup_count = 0;
    bool is_dup = false, is_dup_within = false, is_dup_between = false;
    int32_t j = -1;
    std::vector<uint64_t> nuniq_tiles(ntiles, 0);
    std::vector<uint64_t> ndup_tiles(ntiles, 0);
    std::vector<uint64_t> ndup_tiles_within(ntiles, 0);
    std::vector<uint64_t> ndup_tiles_between(ntiles, 0);
    std::map<uint64_t, uint64_t> dup_hist;
    while( nt5min != UINT64_MAX ) {
        // process the current nt5min
        if ( (nuniq + ndups) % 10000000 == 0 )
            notice("Processing nuniq = %llu, ndups = %llu, ndups_uniq = %llu", nuniq, ndups, ndups_uniq);

        is_dup_within = false;
        if ( imins.size() == 1 ) { // the barcode is probably unique, unless duplicate found in the same tile
            is_dup = false;
            is_dup_between = false;
            dup_count = 1;
            j = imins[0];
            uint64_t next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
            while( next_nt5 == nt5min ) {
                is_dup = true;
                is_dup_within = true;
                ++dup_count;
                next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
            }
            tseqs[j] = next_nt5;
        }
        else { // the barcode is definitely duplicate, across multiple tiles
            is_dup = true;
            is_dup_between = true;
            dup_count = (uint64_t)imins.size();
            for (int32_t i = 0; i < (int32_t)imins.size(); ++i) {
                j = imins[i];
                uint64_t next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
                while( next_nt5 == nt5min ) {
                    is_dup_within = true;
                    ++dup_count;
                    next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
                }
                tseqs[j] = next_nt5;
            }
        }
        // update the statistics and output
        if ( is_dup ) {
           hprintf(wdups, "%s\t%zu\t%llu\n", seqmin.c_str(), imins.size(), dup_count);
        }

        if ( dup_count > 1 )
            ++dup_hist[dup_count];
        if ( is_dup ) {
            ++ndups_uniq;
            ndups += dup_count;
            for (int32_t i = 0; i < (int32_t)imins.size(); ++i) {
                ++ndup_tiles[imins[i]];
            }
            if ( is_dup_within ) {
                ++ndups_within;
                for (int32_t i = 0; i < (int32_t)imins.size(); ++i) {
                    ++ndup_tiles_within[imins[i]];
                }
            }
            if ( is_dup_between ) {
                ++ndups_between;
                for (int32_t i = 0; i < (int32_t)imins.size(); ++i) {
                    ++ndup_tiles_between[imins[i]];
                }
            }
        } else {
            ++nuniq;
            for (int32_t i = 0; i < (int32_t)imins.size(); ++i) {
                ++nuniq_tiles[imins[i]];
            }
        }

        // update the minimum
        nt5min = UINT64_MAX;
        imins.clear();
        for (int32_t i = 0; i < ntiles; ++i)
        {
            if ( nt5min > tseqs[i] ) {
                nt5min = tseqs[i];
                imins.clear();
                imins.push_back(i);
            } else if ( nt5min == tseqs[i] ) {
                imins.push_back(i);
            }
        }
        if ( nt5min != UINT64_MAX ) 
            seqmin.assign(bcdfs[imins[0]]->str_field_at(0)); // assign the minimum sequence
    }
    hts_close(wdups);

    notice("Finished processing nuniq = %llu, ndups = %llu, ndups_uniq = %llu, ndups_uniq_within = %llu, ndups_uniq_between = %llu", nuniq, ndups, ndups_uniq, ndups_within, ndups_between);

    for (int32_t i = 0; i < ntiles; ++i)
    {
        delete bcdfs[i];
    }

    // write the histogram
    sprintf(buf, "%s.hist.tsv", outprefix.c_str());
    htsFile *whist = hts_open(buf, "w");
    if (whist == NULL)
        error("Cannot open %s for writing", buf);
    hprintf(whist, "#dupcount\tnum\n");
    hprintf(whist, "1\t%llu\n", nuniq);
    for(std::map<uint64_t,uint64_t>::iterator it = dup_hist.begin(); it != dup_hist.end(); ++it) {
        hprintf(whist, "%llu\t%llu\n", it->first, it->second);
    }
    hts_close(whist);

    // write the tile counts
    sprintf(buf, "%s.tiles.tsv", outprefix.c_str());
    htsFile *wtile = hts_open(buf, "w");
    if (wtile == NULL)
        error("Cannot open %s for writing", buf);
    hprintf(wtile, "#tile\ttotal\tuniq\tdups_uniq\tdups_uniq_within\tdups_uniq_between\n");
    for (int32_t i = 0; i < df.nrows; ++i)
    {
        hprintf(wtile, "%s\t%llu\t%llu\t%llu\t%llu\t%llu\n", df.get_str_elem(i, "id").c_str(), ntotal_tiles[i], nuniq_tiles[i], ndup_tiles[i], ndup_tiles_within[i], ndup_tiles_between[i]);
    }
    hts_close(wtile);

    notice("Analysis finished");

    return 0;
}
