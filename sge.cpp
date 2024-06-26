#include "sge.h"

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
    //if (!bcd_tr.close())
    //    error("Error in closing the barcode file");
    // if (!ftr_tr.close())
    //     error("Error in closing the feature file");
    // if (!mtx_tr.close())
    //     error("Error in closing the matrix file");
    for (int32_t i = 0; i < (int32_t)ftrs.size(); ++i)
        delete ftrs[i];
    ftrs.clear();
}

bool sge_stream_reader::read_mtx()
{
    // read a line from mtx file
    if (mtx_tr.read_line() == 0)
    { // EOF reached
        notice("EOF reached: cur_line = %llu, cur_iftr = %llu, cur_sbcd.nid = %llu", cur_line, cur_iftr, cur_sbcd.nid);
        return false;
    }
    is_bcd_new = false;
    // parse barcodes until the right barcode was found
    while (cur_sbcd.nid < mtx_tr.uint64_field_at(1))
    {
        if (bcd_tr.read_line() == 0)
        {
            error("EOF reached while finding the barcode %d in the barcode file", mtx_tr.int_field_at(0));
        }
        is_bcd_new = true;
        cur_sbcd.assign_fields(bcd_tr);
    }
    if (cur_sbcd.nid > mtx_tr.uint64_field_at(1))
    {
        error("Cannot find the barcode %llu in the barcode file - cur_sbcd = %s:%d", mtx_tr.uint64_field_at(1), cur_sbcd.strid.c_str(), cur_sbcd.nid);
    }
    if (nfields == 0)
    {
        nfields = mtx_tr.nfields - 2;
    }
    else if (nfields != mtx_tr.nfields - 2)
    {
        error("Inconsistent number of fields in the mtx file - %d vs %d", nfields, mtx_tr.nfields - 2);
    }

    // read the contents of the mtx file
    cur_iftr = mtx_tr.uint64_field_at(0);
    if (cur_cnts.empty())
        cur_cnts.resize(mtx_tr.nfields - 2);
    for (int32_t i = 0; i < mtx_tr.nfields - 2; ++i)
    {
        cur_cnts[i] = mtx_tr.uint64_field_at(i + 2);
    }
    ++cur_line;
    return true;
}

// load the features and fill in the vectors
int32_t sge_stream_reader::load_features()
{
    char buf[65535];
    while (ftr_tr.read_line())
    {
        sge_ftr_t *pftr = new sge_ftr_t(ftr_tr.str_field_at(0),
                                        ftr_tr.str_field_at(1),
                                        ftr_tr.uint64_field_at(2),
                                        ftr_tr.str_field_at(3));
        // ftr.id.assign(ftr_tr.str_field_at(0));
        // ftr.name.assign(ftr_tr.str_field_at(1));
        // ftr.nid = (uint32_t)ftr_tr.int_field_at(2);
        // ftr.cnts.clear();
        // strcpy(buf, ftr_tr.str_field_at(3));
        // const char* pch = buf;
        // while(pch != NULL) {
        //     uint64_t x = (uint64_t)strtoull(pch, NULL, 10);
        //     ftr.cnts.push_back(x);
        //     pch = strchr(pch, ',');
        //     if ( pch != NULL ) ++pch;
        // }
        // ftrs.push_back(ftr);
        ftrs.push_back(pftr);
        if (pftr->nid != (int32_t)ftrs.size())
        {
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
    snprintf(buf, 65535, "%s.tmp", mtxf);
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

    if (cur_sbcd.cnts.empty())
    {
        cur_sbcd.cnts.resize(nfields, 0);
    }
    for (int32_t i = 0; i < nfields; ++i)
    {

        cur_sbcd.cnts[i] = 0;
    }
    return true;
}

bool sge_stream_writer::write_ftr(const char *id, const char *name, uint64_t nid, std::vector<uint64_t> &cnts)
{
    if (cnts.empty())
    {
        cnts.resize(nfields, 0);
    }
    std::string strcnt;
    cat_join_uint64(strcnt, cnts, ",");
    hprintf(wh_ftr, "%s\t%s\t%llu\t%s\n", id, name, nid, strcnt.c_str());
    return true;
}

bool sge_stream_writer::add_mtx(uint64_t iftr, std::vector<uint64_t> &cnts)
{
    if (nfields == 0)
    {
        nfields = (int32_t)cnts.size();
        cur_sbcd.cnts.resize(nfields, 0);
    }
    else if (nfields != (int32_t)cnts.size())
    {
        error("The number of fields in the mtx file is not consistent: %d vs %d", nfields, (int32_t)cnts.size());
    }

    std::string strcnt;
    cat_join_uint64(strcnt, cnts, " ");
    hprintf(wh_tmp, "%llu %llu %s\n", iftr, cur_sbcd.nid, strcnt.c_str());

    // update ftr_cnts
    if (iftr > ftr_cnts.size() + 1)
        ftr_cnts.resize(iftr);
    if (ftr_cnts[iftr - 1].empty())
        ftr_cnts[iftr - 1].resize(nfields, 0);
    if ((ftr_cnts[iftr - 1].size() != nfields) || (cnts.size() != nfields))
        error("The number of fields in the mtx file is not consistent: %d vs %d vs %d", nfields, (int32_t)cnts.size(), (int32_t)ftr_cnts[iftr - 1].size());
    for (int32_t i = 0; i < nfields; ++i)
        ftr_cnts[iftr - 1][i] += cnts[i];

    // update sbcd_cnts
    // assert(cur_sbcd.cnts.size() == nfields);
    for (int32_t i = 0; i < nfields; ++i)
        cur_sbcd.cnts[i] += cnts[i];

    ++nlines;
    return true;
}

// combine header and contents for mtx to a single compressed mtx.gz file
bool sge_stream_writer::flush_mtx()
{
    flush_cur_sbcd();           // write the last barcode
    if (hts_close(wh_tmp) != 0) // close the temporary mtx file
        error("Cannot close the temporary mtx file");
    wh_tmp = NULL;
    if (hts_close(wh_bcd) != 0) // close the barcode file
        error("Cannot close the barcode file");
    wh_bcd = NULL;

    // write a header file
    char buf[65535];
    snprintf(buf, 65535, "%s.hdr", fn_mtx.c_str());
    htsFile *wh_hdr = hts_open(buf, "w"); // write the contents of the matrix file first without the header
    hprintf(wh_hdr, "%%%%MatrixMarket matrix coordinate integer general\n%%\n");
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
    snprintf(buf, 65535, "%s.hdr", fn_mtx.c_str());
    if (remove(buf) != 0)
        error("Cannot remove %s", buf);
    snprintf(buf, 65535, "%s.tmp", fn_mtx.c_str());
    if (remove(buf) != 0)
        error("Cannot remove %s", buf);

    return true;
}

// read a minmax file
bool read_minmax(const char *fn, uint64_t &xmin, uint64_t &xmax, uint64_t &ymin, uint64_t &ymax)
{
    tsv_reader tr(fn);
    if (!tr.read_line())
        error("Cannot read the first line from %s", fn);

    bool xmin_found = false, xmax_found = false, ymin_found = false, ymax_found = false;

    // check the number of fields
    if ( tr.nfields == 2 ) { // two fields - each with xmin, xmax, ymin, ymax
        do {
            if ( tr.nfields != 2 )
                error("The number of fields in the tall-format minmax file is not 2: %d", tr.nfields);
            if (strcmp(tr.str_field_at(0), "xmin") == 0) {
                xmin = tr.uint64_field_at(1);
                xmin_found = true;
            }
            else if (strcmp(tr.str_field_at(0), "xmax") == 0) {
                xmax = tr.uint64_field_at(1);
                xmax_found = true;
            }
            else if (strcmp(tr.str_field_at(0), "ymin") == 0) {
                ymin = tr.uint64_field_at(1);
                ymin_found = true;
            }
            else if (strcmp(tr.str_field_at(0), "ymax") == 0) {
                ymax = tr.uint64_field_at(1);
                ymax_found = true;
            }
            if ( xmin_found && xmax_found && ymin_found && ymax_found )
                break;
        } while ( tr.read_line() );
    }
    else if ( tr.nfields >= 4 ) {
        // check the header line
        int32_t i_xmin = -1, i_xmax = -1, i_ymin = -1, i_ymax = -1;
        for(int32_t i = 0; i < tr.nfields; ++i) {
            if ( strcmp(tr.str_field_at(i), "xmin") == 0 ) i_xmin = i;
            else if ( strcmp(tr.str_field_at(i), "xmax") == 0 ) i_xmax = i;
            else if ( strcmp(tr.str_field_at(i), "ymin") == 0 ) i_ymin = i;
            else if ( strcmp(tr.str_field_at(i), "ymax") == 0 ) i_ymax = i;
        }
        // if all the fields are found
        if ( i_xmin >= 0 && i_xmax >= 0 && i_ymin >= 0 && i_ymax >= 0 ) { // header was found
            if (!tr.read_line()) // read the second line
                error("Cannot read the first second from %s", fn);
            xmin = tr.uint64_field_at(i_xmin);
            xmax = tr.uint64_field_at(i_xmax);
            ymin = tr.uint64_field_at(i_ymin);
            ymax = tr.uint64_field_at(i_ymax);
        }
        else if ( i_xmin < 0 && i_xmax < 0 && i_ymin < 0 && i_ymax < 0 ) { // headerless format
            xmin = tr.uint64_field_at(0);
            xmax = tr.uint64_field_at(1);
            ymin = tr.uint64_field_at(2);
            ymax = tr.uint64_field_at(3);
        }
        else {
            error("Cannot recognize header in the minmax file: %s", fn);
            return false;
        }
    }
    else {
        error("The number of fields in the minmax file is not recognized: %d", tr.nfields);
        return false;
    }
    return true;
}