#include "tiles.h"

///////////////////////////////////////////////////////////
// METHODS FOR tile_writer
///////////////////////////////////////////////////////////
tile_writer::tile_writer(const char *_rootdir, const char *_filename, bool _gzip) : gzip(_gzip)
{
    rootdir.assign(_rootdir);
    if (rootdir[rootdir.size() - 1] != '/')
        rootdir += "/";
    filename.assign(_filename);
    all_filename = rootdir + filename;
    all_fh = hts_open(all_filename.c_str(), gzip ? "wz" : "w");
}

void tile_writer::close()
{
    if (all_fh != NULL)
    {
        hts_close(all_fh);
        all_fh = NULL;
    }
    for (std::map<int32_t, htsFile *>::iterator it = lane_fhs.begin(); it != lane_fhs.end(); ++it)
    {
        hts_close(it->second);
    }
    lane_fhs.clear();

    for (std::map<int32_t, std::map<int32_t, htsFile *>>::iterator it = tile_fhs.begin(); it != tile_fhs.end(); ++it)
    {
        for (std::map<int32_t, htsFile *>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            hts_close(it2->second);
        }
    }
    tile_fhs.clear();
}

bool tile_writer::write_lane(int32_t lane, const std::string &s)
{
    std::map<int32_t, htsFile *>::iterator it = lane_fhs.find(lane);
    bool ret = false;
    htsFile *wf = NULL;
    if (it == lane_fhs.end())
    {
        std::string name;
        catprintf(name, "%s%d", rootdir.c_str(), lane);
        ret = makePath(name); // make directory if it does not exist
        catprintf(name, "/%s", filename.c_str());
        wf = hts_open(name.c_str(), gzip ? "wz" : "w");
        lane_fhs[lane] = wf;
        lane_filenames[lane] = name;
    }
    else
    {
        wf = it->second;
    }
    hprint_str(wf, s);
    return ret;
}

bool tile_writer::write_tile(int32_t lane, int32_t tile, const std::string &s)
{
    std::map<int32_t, htsFile *>::iterator it = tile_fhs[lane].find(tile);
    bool ret = false;
    htsFile *wf = NULL;
    if (it == tile_fhs[lane].end())
    {
        std::string name;
        catprintf(name, "%s%d/%d", rootdir.c_str(), lane, tile);
        ret = makePath(name); // make directory if it does not exist
        catprintf(name, "/%s", filename.c_str());
        wf = hts_open(name.c_str(), gzip ? "wz" : "w");
        tile_fhs[lane][tile] = wf;
        tile_filenames[lane][tile] = name;
    }
    else
    {
        wf = it->second;
    }
    hprint_str(wf, s);
    return ret;
}

/////////////////////////////////////////////////////////////
// METHODS FOR tile_counter
/////////////////////////////////////////////////////////////

bool tile_counter::add_count(int32_t lane, int32_t tile, int32_t icol, int32_t cnt)
{
    all_cnts[icol] += cnt;

    std::vector<uint64_t> &lcnts = lane_cnts[lane];
    if (lcnts.empty())
    {
        lcnts.resize(ncols, 0);
        lcnts[icol] = cnt;
    }
    else
    {
        lcnts[icol] += cnt;
    }

    std::vector<uint64_t> &tcnts = tile_cnts[lane][tile];
    if (tcnts.empty())
    {
        tcnts.resize(ncols, 0);
        tcnts[icol] = cnt;
        return true;
    }
    else
    {
        tcnts[icol] += cnt;
        return false;
    }
}

bool tile_counter::add_counts(int32_t lane, int32_t tile, std::vector<int32_t> &v)
{
    for (int32_t i = 0; i < ncols; ++i)
        all_cnts[i] += v[i];

    std::vector<uint64_t> &lcnts = lane_cnts[lane];
    if (lcnts.empty())
    {
        lcnts.resize(ncols, 0);
        for (int32_t i = 0; i < ncols; ++i)
            lcnts[i] = v[i];
    }
    else
    {
        for (int32_t i = 0; i < ncols; ++i)
            lcnts[i] += v[i];
    }

    std::vector<uint64_t> &tcnts = tile_cnts[lane][tile];
    if (tcnts.empty())
    {
        tcnts.resize(ncols, 0);
        for (int32_t i = 0; i < ncols; ++i)
            tcnts[i] = v[i];
        return true;
    }
    else
    {
        for (int32_t i = 0; i < ncols; ++i)
            tcnts[i] += v[i];
        return false;
    }
}

/////////////////////////////////////////////////////////////
// METHODS FOR sbcd_sync_reader
/////////////////////////////////////////////////////////////
tsv_reader *sbcd_sync_reader::move_to_ibcd(uint64_t ibcd)
{
    // find the right barcode first
    while (bcd_tr.nfields > 0 && bcd_tr.nlines < ibcd)
    {
        bcd_tr.read_line();
    }
    bcd = bcd_tr.str_field_at(0);
    lbcd = strlen(bcd);

    // scan the matching barcodes from the sbcd files
    int32_t nsbcds = (int32_t)sbcd_trs.size();
    int32_t cmp;
    idx_match = -1;
    for (int32_t i = 0; i < nsbcds; ++i)
    {
        if (sbcd_trs[i]->nfields > 0)
        { // consider only when EOF is not yet reached
            while ((cmp = strncmp(bcd, sbcd_trs[i]->str_field_at(0), lbcd)) > 0)
            {
                if (sbcd_trs[i]->read_line() == 0)
                {
                    cmp = -1;
                    break;
                }
            }
            if (cmp == 0)
            {
                // match found
                idx_match = i;
                return sbcd_trs[idx_match];
                // return true;
            }
        }
    }
    return NULL;
    // return false; // match not found
}

/////////////////////////////////////////////////////////////
// Global method
/////////////////////////////////////////////////////////////
// once a single barcode finishes reading, flush the barcode out to all, lane, tile
void write_sbcd(sbcd_sync_reader &ssr, uint64_t cur_ibcd, std::vector<int32_t> &cur_bcd_cnts, tile_writer &bcd_tw, tile_counter &sbcds_counter, int32_t n_mtx)
{
    // print the current barcode
    tsv_reader *p = ssr.move_to_ibcd(cur_ibcd);
    int32_t lane = p == NULL ? 0 : p->int_field_at(1);
    int32_t tile = p == NULL ? 0 : p->int_field_at(2);
    int32_t x = p == NULL ? 0 : p->int_field_at(3);
    int32_t y = p == NULL ? 0 : p->int_field_at(4);

    std::string outstr;

    // write all
    catprintf(outstr, "%s\t%llu\t%llu\t%d\t%d\t%d\t%d\t", ssr.bcd, sbcds_counter.all_cnts[0] + 1, cur_ibcd, lane, tile, x, y);
    cat_join_int32(outstr, cur_bcd_cnts, ",");
    outstr += "\n";
    bcd_tw.write_all(outstr);

    outstr.clear();

    // write lane
    catprintf(outstr, "%s\t%llu\t%llu\t%d\t%d\t%d\t%d\t", ssr.bcd, sbcds_counter.get_lane_count(lane, 0) + 1, cur_ibcd, lane, tile, x, y);
    cat_join_int32(outstr, cur_bcd_cnts, ",");
    outstr += "\n";
    bcd_tw.write_lane(lane, outstr);

    outstr.clear();

    // write tile
    catprintf(outstr, "%s\t%llu\t%llu\t%d\t%d\t%d\t%d\t", ssr.bcd, sbcds_counter.get_tile_count(lane, tile, 0) + 1, cur_ibcd, lane, tile, x, y);
    cat_join_int32(outstr, cur_bcd_cnts, ",");
    outstr += "\n";
    bcd_tw.write_tile(lane, tile, outstr);

    sbcds_counter.add_count(lane, tile, 0, 1); // add # barcodes
}


void open_tiles(std::vector<std::string> &tile_paths, std::vector<tsv_reader*>& bcdfs) {
    // clear up the existing files if exists
    if (bcdfs.size() > 0)
    {
        for (int32_t i = 0; i < (int32_t)bcdfs.size(); ++i)
        {
            delete bcdfs[i];
        }
        bcdfs.clear();
    }

    for(int32_t i=0; i < (int32_t)tile_paths.size(); ++i) {
        bcdfs.push_back(new tsv_reader(tile_paths[i].c_str()));
        if (bcdfs.back()->read_line() == 0)
        {
            error("ERROR: Observed an empty barcode file %s", tile_paths[i].c_str());
        }
    }
}

// open all tiles
void open_tiles(dataframe_t &df, std::vector<std::string> &tiles, std::vector<tsv_reader *> &bcdfs)
{
    // read manifest files
    if (bcdfs.size() > 0)
    {
        for (int32_t i = 0; i < (int32_t)bcdfs.size(); ++i)
        {
            delete bcdfs[i];
        }
        bcdfs.clear();
    }

    tiles = df.get_column("id");
    int32_t icol = df.get_colidx("fullpath");
    for (int32_t i = 0; i < df.nrows; ++i)
    {
        bcdfs.push_back(new tsv_reader(df.get_str_elem(i, icol).c_str()));
        if (bcdfs.back()->read_line() == 0)
        {
            error("ERROR: Observed an empty barcode file %s", df.get_str_elem(i, icol).c_str());
        }
    }
}

std::pair<uint64_t, uint64_t> count_matches(std::vector<uint64_t> &bseqs, dataframe_t &df, std::vector<uint64_t> &dcounts, int32_t match_len, htsFile *wmatch)
{
    std::vector<std::string> tiles;
    std::vector<tsv_reader *> bcdfs;

    open_tiles(df, tiles, bcdfs);

    int32_t ntiles = (int32_t)tiles.size();
    if (dcounts.empty())
    {
        dcounts.resize(ntiles, 0);
    }
    int32_t len = strlen(bcdfs[0]->str_field_at(0));
    if (len < match_len)
        error("HDMI length %d does not match to the parameters %d", len, match_len);

    std::vector<uint64_t> tseqs(ntiles);
    for (int32_t i = 0; i < ntiles; ++i)
    {
        tseqs[i] = seq2nt5(bcdfs[i]->str_field_at(0), match_len);
    }

    // sort the batch of sequences
    uint64_t batch_size = (uint64_t)bseqs.size();
    notice("Started sorting of %llu records", batch_size);
    std::sort(bseqs.begin(), bseqs.end());
    notice("Finished sorting of %llu records", batch_size);
    uint64_t nseqs = (uint64_t)bseqs.size();

    uint64_t nmiss = 0;
    uint64_t ndups = 0;
    bool has_match, is_dup;
    int32_t cmp;
    // count the sequences that matches
    for (uint64_t i = 0; i < nseqs; ++i)
    {
        if (i % (batch_size / 20) == 0)

            notice("Processing %d records, nmiss = %llu, ndups = %llu, bseqs[i]=%032llu, tseqs[0]=%032llu", i, nmiss, ndups, bseqs[i], tseqs[0]);

        has_match = false;
        is_dup = false;
        uint64_t s = bseqs[i];
        for (int32_t j = 0; j < ntiles; ++j)
        {
            cmp = s < tseqs[j] ? -1 : (s == tseqs[j] ? 0 : 1);
            while (cmp > 0)
            {
                if (bcdfs[j]->read_line() == 0)
                    tseqs[j] = UINT64_MAX;
                else
                    tseqs[j] = seq2nt5(bcdfs[j]->str_field_at(0), match_len);
                cmp = s < tseqs[j] ? -1 : (s == tseqs[j] ? 0 : 1);
            }
            if (cmp == 0)
            {
                has_match = true;
                hprintf(wmatch, "%s", bcdfs[j]->str_field_at(0));
                for (int32_t k = 1; k < bcdfs[j]->nfields; ++k)
                    hprintf(wmatch, "\t%s", bcdfs[j]->str_field_at(k));
                hprintf(wmatch, "\n");
                ++dcounts[j];
                if ((i > 0) && (bseqs[i] == bseqs[i - 1]))
                {
                    is_dup = true;
                }
                /*
                else {
                  ++ucounts[j];
                }
                */
            }
        }
        if (is_dup)
        {
            ++ndups;
        }
        else if (!has_match)
        {
            ++nmiss;
        }
    }
    notice("Finished processing a batch of %d records, nmiss = %llu, ndups = %llu", nseqs, nmiss, ndups);

    for (int32_t i = 0; i < ntiles; ++i)
    {
        delete bcdfs[i];
    }

    // return nmiss + ndups;
    return std::make_pair(nmiss, ndups);
}

uint64_t read_bcdf(tsv_reader* bcdf, int32_t match_len, uint64_t& cnt) {
    if ( bcdf->read_line() == 0 ) return UINT64_MAX;
    else {
        ++cnt;
        return seq2nt5(bcdf->str_field_at(0), match_len);
    }
}

std::pair<uint64_t, uint64_t> count_matches_skip_dups(std::vector<uint64_t> &bseqs, dataframe_t &df, std::vector<uint64_t> &dcounts, int32_t match_len, htsFile *wmatch)
{
    std::vector<std::string> tiles;
    std::vector<tsv_reader *> bcdfs;

    open_tiles(df, tiles, bcdfs);

    int32_t ntiles = (int32_t)tiles.size();
    if (dcounts.empty())
    {
        dcounts.resize(ntiles, 0);
    }
    int32_t len = strlen(bcdfs[0]->str_field_at(0));
    if (len < match_len)
        error("HDMI length %d does not match to the parameters %d", len, match_len);

    // sort the batch of sequences
    uint64_t batch_size = (uint64_t)bseqs.size();
    notice("Started sorting of %llu records", batch_size);
    std::sort(bseqs.begin(), bseqs.end());
    notice("Finished sorting of %llu records", batch_size);
    uint64_t nseqs = (uint64_t)bseqs.size();

    size_t cur_bseq = 0; // index for bseq

    // read the first entry in each tile, and store the minimum values and locations
    std::vector<uint64_t> tseqs(ntiles);
    std::vector<int32_t> imins;
    uint64_t nt5min = UINT64_MAX;
    std::string seqmin;
    std::vector<uint64_t> valmins;
    std::vector<uint64_t> ntotal_tiles(ntiles, 0);
    for (int32_t i = 0; i < ntiles; ++i)
    {
        tseqs[i] = read_bcdf(bcdfs[i], match_len, ntotal_tiles[i]);
        if (nt5min > tseqs[i])
        {
            nt5min = tseqs[i];
            imins.clear();
            imins.push_back(i);
        }
        else if (nt5min == tseqs[i])
        {
            imins.push_back(i);
        }
    }
    seqmin.assign(bcdfs[imins[0]]->str_field_at(0)); // assign the minimum sequence
    valmins.resize(bcdfs[imins[0]]->nfields);
    for (int32_t i = 1; i < bcdfs[imins[0]]->nfields; ++i)
    {
        valmins[i] = bcdfs[imins[0]]->uint64_field_at(i);
    }

    // keeping track of different types of duplicates:
    // 1. 2nd-seq duplicate barcodes : multiple 2nd-seq barcodes observed (not really necessary)
    // 2. 1st-seq duplicate barcodes : multiple 1st-seq barcodes observed
    uint64_t nmatch_uniq_1st = 0, nmatch_dups_1st = 0, nskip_uniq_1st = 0, nskip_dups_1st = 0;
    uint64_t nmatch_uniq_2nd = 0, nmatch_dups_2nd = 0, nskip_2nd = 0;
    bool is_dup = false;
    int32_t j = -1;
    std::vector<uint64_t> nmatch_uniq_tiles(ntiles, 0); // number of unique barcodes in each tile
    std::vector<uint64_t> nmatch_dup_tiles(ntiles, 0);  // number of duplicate barcodes in each tile
    while (nt5min != UINT64_MAX) // check the minimum value is still valid
    {
        // process the current nt5min
        if ((nmatch_uniq_1st + nmatch_dups_1st + nskip_uniq_1st + nskip_dups_1st) % 10000000 == 0)
            notice("Processing nmatch_uniq_1st = %llu, nmatch_dups_1st = %llu, nskip_uniq_1st = %llu, nskip_dups_1st = %llu, nmatch_uniq_2nd = %llu, nmatch_dups_2nd = %llu, nskip_2nd = %llu", 
                nmatch_uniq_1st, nmatch_dups_1st, nskip_uniq_1st, nskip_dups_1st, nmatch_uniq_2nd, nmatch_dups_2nd, nskip_2nd);

        if (imins.size() == 1)
        { // the barcode is probably unique, unless duplicate found in the same tile
            is_dup = false;
            j = imins[0];
            uint64_t next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
            while (next_nt5 == nt5min)
            {
                is_dup = true;
                next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
            }
            tseqs[j] = next_nt5;
        }
        else
        { // the barcode is definitely duplicate, across multiple tiles
            is_dup = true;
            for (int32_t i = 0; i < (int32_t)imins.size(); ++i)
            {
                j = imins[i];
                uint64_t next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
                while (next_nt5 == nt5min)
                {
                    next_nt5 = read_bcdf(bcdfs[j], match_len, ntotal_tiles[j]);
                }
                tseqs[j] = next_nt5;
            }
        }

        // count non-matches first
        while( bseqs[cur_bseq] < nt5min ) { 
            ++nskip_2nd;
            ++cur_bseq;
        }
        if ( bseqs[cur_bseq] == nt5min ) { // there is a match
            // write out the matching barcodes if it is not a duplicate
            if ( is_dup ) {
                ++nmatch_dups_1st;
            }
            else {  // only write unique matches here, only once per barcode
                ++nmatch_uniq_1st;
                j = imins[0];
                // write out the barcodes
                if ( bcdfs[j]->nfields > 1) {
                    hprintf(wmatch, "%s", seqmin.c_str());
                    for (int32_t k = 1; k < (int32_t)valmins.size(); ++k)
                        hprintf(wmatch, "\t%llu", valmins[k]);
                    // hprintf(wmatch, "%s", bcdfs[j]->str_field_at(0));
                    // for (int32_t k = 1; k < bcdfs[j]->nfields; ++k)
                    //     hprintf(wmatch, "\t%s", bcdfs[j]->str_field_at(k));
                    hprintf(wmatch, "\n");
                }
            }

            // count the 2nd-seq duplicates
            bool is_dup_2nd = false;
            while( bseqs[cur_bseq] == nt5min ) {
                if (is_dup)
                {
                    ++nmatch_dups_2nd;
                    for (int32_t i = 0; i < (int32_t)imins.size(); ++i)
                    {
                        ++nmatch_dup_tiles[imins[i]];
                    }
                }
                else
                {
                    ++nmatch_uniq_2nd;
                    for (int32_t i = 0; i < (int32_t)imins.size(); ++i)
                    {
                        ++nmatch_uniq_tiles[imins[i]];
                    }
                }
                ++cur_bseq;
                if ( bseqs[cur_bseq] == nt5min )
                    is_dup_2nd = true;
            }
        }
        else { // there is no match
            if ( is_dup ) {
                ++nskip_dups_1st;
            }
            else {
                ++nskip_uniq_1st;
            }
        }

        // update the minimum
        nt5min = UINT64_MAX;
        imins.clear();
        for (int32_t i = 0; i < ntiles; ++i)
        {
            if (nt5min > tseqs[i])
            {
                nt5min = tseqs[i];
                imins.clear();
                imins.push_back(i);
            }
            else if (nt5min == tseqs[i])
            {
                imins.push_back(i);
            }
        }
        if (nt5min != UINT64_MAX) {
            seqmin.assign(bcdfs[imins[0]]->str_field_at(0)); // assign the minimum sequence
            //valmins.resize(bcdfs[imins[0]]->nfields);
            for (int32_t i = 1; i < bcdfs[imins[0]]->nfields; ++i)
            {
                valmins[i] = bcdfs[imins[0]]->uint64_field_at(i);
            }
        }
    }

    for (int32_t i = 0; i < ntiles; ++i)
    {
        delete bcdfs[i];
        dcounts[i] += nmatch_uniq_tiles[i];
    }

    // return nmiss + ndups;
    return std::make_pair(nskip_2nd + nmatch_dups_2nd, nmatch_dups_2nd);
}
