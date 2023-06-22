#include "sge.h"

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
// METHODS FOR sge_stream_reader
/////////////////////////////////////////////////////////////
void sge_stream_reader::open(const char *bcdf, const char *ftrf, const char *mtxf)
{
    bcd_tr.open(bcdf);
    ftr_tr.open(ftrf);
    mtx_tr.delimiter = ' ';
    mtx_tr.open(mtxf);

    // read the header of the matrix
    mtx_tr.read_line();
    while (mtx_tr.str_field_at(0)[0] == '%')
    {
        mtx_tr.read_line();
    }
    nftrs = mtx_tr.uint64_field_at(0);
    nbcds = mtx_tr.uint64_field_at(1);
    nlines = mtx_tr.uint64_field_at(2);

    cur_line = cur_iftr = 0;
    cur_sbcd.nid = 0;
    is_bcd_new = false;
    nfields = 0;
}

void sge_stream_reader::close()
{
    bcd_tr.close();
    ftr_tr.close();
    mtx_tr.close();
}

bool sge_stream_reader::read_mtx()
{
    // read a line from mtx file
    if (mtx_tr.read_line() == 0)
    { // EOF reached
        return false;
    }
    is_bcd_new = false;
    // parse barcodes until the right barcode was found
    while (cur_sbcd.nid < mtx_tr.int_field_at(1))
    {
        if (bcd_tr.read_line() == 0)
        {
            error("EOF reached while finding the barcode %lu in the barcode file", mtx_tr.int_field_at(0));
        }
        is_bcd_new = true;
        cur_sbcd.assign_fields(bcd_tr);
    }
    if (cur_sbcd.nid > mtx_tr.int_field_at(1))
    {
        error("Cannot find the barcode %lu in the barcode file - cur_sbcd = %s:%d", mtx_tr.int_field_at(1), cur_sbcd.strid.c_str(), cur_sbcd.nid);
    }
    if ( nfields == 0 ) {
        nfields = mtx_tr.nfields - 2;
    }
    else if ( nfields != mtx_tr.nfields - 2 ) {
        error("Inconsistent number of fields in the mtx file - %d vs %d", nfields, mtx_tr.nfields - 2);
    }

    // read the contents of the mtx file
    cur_iftr = mtx_tr.int_field_at(0);
    cur_cnts.resize(mtx_tr.nfields - 2);
    for (int32_t i = 0; i < cur_cnts.size(); ++i)
    {
        cur_cnts[i] = mtx_tr.uint64_field_at(i + 2);
    }
    ++cur_line;
    return true;
}

// load the features and fill in the vectors
int32_t sge_stream_reader::load_features()
{
    while (ftr_tr.read_line())
    {
        sge_ftr_t ftr;
        ftr.id.assign(ftr_tr.str_field_at(0));
        ftr.name.assign(ftr_tr.str_field_at(1));
        ftr.nid = ftr_tr.int_field_at(2);
        std::vector<std::string> v;
        split(v, ",", ftr_tr.str_field_at(3));
        ftr.cnts.resize(v.size());
        for (int32_t i = 0; i < (int32_t)v.size(); ++i)
        {
            ftr.cnts[i] = atoi(v[i].c_str());
        }
        ftrs.push_back(ftr);
        if ( ftr.nid != (int32_t)ftrs.size() ) {
            error("Feature nid should be 1-based sequential number");
        }
    }
    return (int32_t)ftrs.size();
}

/////////////////////////////////////////////////////////////
// METHODS FOR sge_stream_writer
/////////////////////////////////////////////////////////////
void sge_stream_writer::open(const char *bcdf, const char *ftrf, const char *mtxf)
{
    char buf[65535];
    fn_mtx.assign(mtxf);
    sprintf(buf, "%s.tmp", mtxf);
    wh_bcd = hts_open(bcdf, "wz");
    wh_ftr = hts_open(ftrf, "wz");
    wh_tmp = hts_open(buf, "w"); // write the contents of the matrix file first without the header
    if (wh_bcd == NULL || wh_ftr == NULL || wh_tmp == NULL)
    {
        error("Cannot open the output files");
    }
    cur_sbcd.nid = cur_sbcd.gid = nfields = nlines = 0;
}

void sge_stream_writer::close()
{
    if (wh_tmp != NULL)
    {
        if (!flush_mtx())
        {
            error("Cannot flush the mtx file");
        }
    }
    // if (wh_bcd != NULL)
    // {
    //     if (!hts_close(wh_bcd))
    //     {
    //         error("Cannot close the barcode file");
    //     }
    //     wh_bcd = NULL;
    // }
    if (wh_ftr != NULL)
    {
        if (hts_close(wh_ftr) != 0)
        {
            error("Cannot close the feature file");
        }
        wh_ftr = NULL;
    }
}

bool sge_stream_writer::flush_cur_sbcd()
{
    std::string strcnt;
    cat_join_uint64(strcnt, cur_sbcd.cnts, ",");
    hprintf(wh_bcd, "%s\t%llu\t%llu\t%lu\t%lu\t%llu\t%llu\t%s\n", 
            cur_sbcd.strid.c_str(), cur_sbcd.nid, cur_sbcd.gid, cur_sbcd.lane, cur_sbcd.tile, cur_sbcd.px, cur_sbcd.py, strcnt.c_str());
    // hprintf(wh_bcd, "%s\t", cur_sbcd.strid.c_str()); 
    // hprintf(wh_bcd, "%llu\t", cur_sbcd.nid);  
    // hprintf(wh_bcd, "%llu\t", cur_sbcd.gid);  
    // hprintf(wh_bcd, "%lu\t", cur_sbcd.lane);  
    // hprintf(wh_bcd, "%lu\t", cur_sbcd.tile);  
    // hprintf(wh_bcd, "%llu\t", cur_sbcd.px);
    // hprintf(wh_bcd, "%llu\t", cur_sbcd.py);
    // hprintf(wh_bcd, "%s\n", strcnt.c_str());      
    return true;
}

// add a spatial barcode
bool sge_stream_writer::add_sbcd(const char *strid, uint64_t old_nid, uint32_t lane, uint32_t tile, uint64_t x, uint64_t y)
{
    if (cur_sbcd.nid > 0)
    {
        flush_cur_sbcd();
    }
    cur_sbcd.strid.assign(strid);
    cur_sbcd.nid++;
    cur_sbcd.gid = old_nid;
    cur_sbcd.lane = lane;
    cur_sbcd.tile = tile;
    cur_sbcd.px = x;
    cur_sbcd.py = y;

    if ( cur_sbcd.cnts.empty() ) {
        cur_sbcd.cnts.resize(nfields, 0);
    }
    for(int32_t i=0; i < nfields; ++i) {
     
        cur_sbcd.cnts[i] = 0;
    }   
    return true;
}

bool sge_stream_writer::write_ftr(const char *id, const char *name, uint64_t nid, std::vector<uint64_t>& cnts)
{
    if ( cnts.empty() ) {
        cnts.resize(nfields, 0);
    }
    std::string strcnt;
    cat_join_uint64(strcnt, cnts, ",");
    hprintf(wh_ftr, "%s\t%s\t%llu\t%s\n", id, name, nid, strcnt.c_str());
    return true;
}

bool sge_stream_writer::add_mtx(uint64_t iftr, std::vector<uint64_t> &cnts)
{
    if ( nfields == 0 ) {
        nfields = (int32_t)cnts.size();
        cur_sbcd.cnts.resize(nfields, 0);
    }
    else if ( nfields != (int32_t)cnts.size() ) {
        error("The number of fields in the mtx file is not consistent: %d vs %d", nfields, (int32_t)cnts.size());
    }

    std::string strcnt;
    cat_join_uint64(strcnt, cnts, " ");
    hprintf(wh_tmp, "%llu %llu %s\n", iftr, cur_sbcd.nid, strcnt.c_str());

    // update ftr_cnts
    if ( iftr > ftr_cnts.size() )
        ftr_cnts.resize(iftr+1);
    ftr_cnts[iftr].resize(nfields);
    for (int32_t i = 0; i < nfields; ++i )
        ftr_cnts[iftr][i] += cnts[i];

    // update sbcd_cnts
    for (int32_t i = 0; i < nfields; ++i )
        cur_sbcd.cnts[i] += cnts[i];
    
    ++nlines;
    return true;
}

// combine header and contents for mtx to a single compressed mtx.gz file
bool sge_stream_writer::flush_mtx()
{
    flush_cur_sbcd(); // write the last barcode
    if ( hts_close(wh_tmp) != 0 ) // close the temporary mtx file
        error("Cannot close the temporary mtx file");
    wh_tmp = NULL;
    if ( hts_close(wh_bcd) != 0 ) // close the barcode file
        error("Cannot close the barcode file");
    wh_bcd = NULL;

    // write a header file
    char buf[65535];
    sprintf(buf, "%s.hdr", fn_mtx.c_str());
    htsFile *wh_hdr = hts_open(buf, "w"); // write the contents of the matrix file first without the header
    hprintf(wh_hdr, "%%MatrixMarket matrix coordinate integer general\n%\n");
    hprintf(wh_hdr, "%d %llu %llu\n", (int32_t)ftr_cnts.size(), cur_sbcd.nid, nlines);
    hts_close(wh_hdr);

    // merge the two files, requires gzip
    notice("Generating the merged %s file", fn_mtx.c_str());
    std::string cmd;
    catprintf(cmd, "cat %s.hdr %s.tmp | gzip -c > %s", fn_mtx.c_str(), fn_mtx.c_str(), fn_mtx.c_str());
    int32_t ret = system(cmd.c_str());
    if ((ret == -1) || (WEXITSTATUS(ret) == 127))
    {
        error("Error in running %s", cmd.c_str());
    }

    // remove the header and tmp file
    sprintf(buf, "%s.hdr", fn_mtx.c_str());
    if (remove(buf) != 0)
        error("Cannot remove %s", buf);
    sprintf(buf, "%s.tmp", fn_mtx.c_str());
    if (remove(buf) != 0)
        error("Cannot remove %s", buf);

    return true;
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

// read a minmax file
bool read_minmax(const char *fn, uint64_t &xmin, uint64_t &xmax, uint64_t &ymin, uint64_t &ymax) {
    tsv_reader tr(fn);
    // if ( !tr.is_open() )
    //     error("Cannot open %s", fn);
    if ( !tr.read_line() )
        error("Cannot read the first line from %s", fn);
    if ( tr.nfields != 4 )
        error("The number of fields in the minmax file is not 4: %d", tr.nfields);
    
    if ( strcmp(tr.str_field_at(0),"xmin") == 0 ) { // header exists, so read the second line
        if ( !tr.read_line() )
            error("Cannot read the first line from %s", fn);
    }
    xmin = tr.uint64_field_at(0);
    xmax = tr.uint64_field_at(1);
    ymin = tr.uint64_field_at(2);
    ymax = tr.uint64_field_at(3);
    return true;
}