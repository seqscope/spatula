#include "sge2.h"

/////////////////////////////////////////////////////////////
// METHODS FOR sge2_stream_reader
/////////////////////////////////////////////////////////////
void sge2_stream_reader::open(const char *bcdf, const char *ftrf, const char *mtxf)
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

void sge2_stream_reader::close()
{
    bcd_tr.close();
    ftr_tr.close();
    mtx_tr.close();
    for (int32_t i = 0; i < (int32_t)ftrs.size(); ++i)
        delete ftrs[i];
    ftrs.clear();
}

bool sge2_stream_reader::read_mtx()
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

        // parse the header if exists
        // if ( cur_sbcd.nid == 0 && bcd_tr.str_field_at(0)[0] == '#' ) {
        //     icol_bcd_strid = icol_bcd_px = icol_bcd_py = -1;
        //     for(int32_t i = 0; i < bcd_tr.nfields; ++i) {
        //         if ( strcmp(bcd_tr.str_field_at(i), colname_bcd_strid.c_str()) == 0 ) icol_bcd_strid = i;
        //         else if ( strcmp(bcd_tr.str_field_at(i), colname_bcd_px.c_str()) == 0 ) icol_bcd_px = i;
        //         else if ( strcmp(bcd_tr.str_field_at(i), colname_bcd_py.c_str()) == 0 ) icol_bcd_py = i;
        //     }

        //     if ( icol_bcd_strid < 0 ) {
        //         error("Cannot find the columns %s, in the barcode file", colname_bcd_strid.c_str());
        //     }

        //     if ( bcd_pos_map.empty() && ( icol_bcd_px < 0 || icol_bcd_py < 0 ) ) {
        //         error("Cannot find the columns %s, and %s in the barcode file without pos file", colname_bcd_px.c_str(), colname_bcd_py.c_str());
        //     }
        // }
        // else 
        if ( cur_sbcd.nid == 0 ) { // no header, use default indices
            if ( bcd_pos_map.empty() ) {
                if ( icol_bcd_strid < 0 || icol_bcd_px < 0 || icol_bcd_py < 0 ) {
                    error("Indices for barcode, X, Y are not set properly - %d, %d, %d", icol_bcd_strid, icol_bcd_px, icol_bcd_py);
                }
                if ( icol_bcd_strid >= bcd_tr.nfields || icol_bcd_px >= bcd_tr.nfields || icol_bcd_py >= bcd_tr.nfields ) {
                    error("The number of fields %d in the barcode file is less than the indices %d, %d, and %d", bcd_tr.nfields, icol_bcd_strid, icol_bcd_px, icol_bcd_py);
                }
            }
            else {
                if ( icol_bcd_strid < 0 ) {
                    error("Negative icol_bcd_strid = %d", icol_bcd_strid);
                }
                if ( icol_bcd_strid >= bcd_tr.nfields ) {
                    error("The number of fields in the barcode file is less than the indices %d", icol_bcd_strid);
                }
                icol_bcd_px = icol_bcd_py = -1; // ignore the position columns
            }
        }
        const char* strid = bcd_tr.str_field_at(icol_bcd_strid);
        if ( bcd_pos_map.empty() ) {
            double px = bcd_tr.double_field_at(icol_bcd_px);
            double py = bcd_tr.double_field_at(icol_bcd_py);
            cur_sbcd.update_values(strid, px, py);
        }
        else {
            auto it = bcd_pos_map.find(strid);
            if ( it == bcd_pos_map.end() ) {
                error("Cannot find the barcode %s in the position file", strid);
            }
            cur_sbcd.update_values(strid, it->second.first, it->second.second);
        }
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
int32_t sge2_stream_reader::load_features()
{
    char buf[65535];
    while (ftr_tr.read_line())
    {
        // header line exists
        // if ( ftrs.empty() && ftr_tr.str_field_at(0)[0] == '#' ) {
        //     icol_ftr_id = icol_ftr_name = -1;
        //     for(int32_t i = 0; i < ftr_tr.nfields; ++i) {
        //         if ( strcmp(ftr_tr.str_field_at(i), colname_ftr_id.c_str()) == 0 ) icol_ftr_id = i;
        //         else if ( strcmp(ftr_tr.str_field_at(i), colname_ftr_name.c_str()) == 0 ) icol_ftr_name = i;
        //     }
        //     if ( icol_ftr_id < 0 || icol_ftr_name < 0 ) {
        //         error("Cannot find the columns %s and %s in the feature file", colname_ftr_id.c_str(), colname_ftr_name.c_str());
        //     }
        // }
        // else {

        if ( icol_ftr_id < 0 || icol_ftr_name < 0 ) {
            error("icol_ftr_id = %d, icol_ftr_name = %d must be non-negative", icol_ftr_id, icol_ftr_name);
        }
        if ( icol_ftr_id >= ftr_tr.nfields || icol_ftr_name >= ftr_tr.nfields ) {
            error("icol_ftr_id = %d, icol_ftr_name = %d must be less than nfields = %d", icol_ftr_id, icol_ftr_name, ftr_tr.nfields);
        }

        sge2_ftr_t *pftr = new sge2_ftr_t(ftr_tr.str_field_at(icol_ftr_id), ftr_tr.str_field_at(icol_ftr_name));
        pftr->nid = (int32_t)ftrs.size() + 1;
        ftrs.push_back(pftr);
        // if (pftr->nid != (int32_t)ftrs.size())
        // {
        //     error("Feature nid should be 1-based sequential number");
        // }
        //}
    }
    return (int32_t)ftrs.size();
}

uint64_t sge2_stream_reader::load_position_file(const char* filename, const char* colname_strid, const char* colname_px, const char* colname_py, char sep) {
    tsv_reader pos_tr;
    pos_tr.delimiter = (int32_t)sep;
    pos_tr.open(filename);
    
    // read the header to identify the columns
    int32_t icol_strid = -1;
    int32_t icol_px = -1;
    int32_t icol_py = -1;

    if (pos_tr.read_line()) {
        for(int32_t i = 0; i < pos_tr.nfields; ++i) {
            if ( strcmp(pos_tr.str_field_at(i), colname_strid) == 0 ) icol_strid = i;
            else if ( strcmp(pos_tr.str_field_at(i), colname_px) == 0 ) icol_px = i;
            else if ( strcmp(pos_tr.str_field_at(i), colname_py) == 0 ) icol_py = i;
        }
    }
    else {
        error("Cannot read the first line from %s", filename);
    
    }

    if ( icol_strid < 0 || icol_px < 0 || icol_py < 0 ) {
        error("Cannot find the columns %s, %s, and %s in the position file %s", colname_strid, colname_px, colname_py, filename);
    }

    return load_position_file_helper(pos_tr, icol_strid, icol_px, icol_py, sep);
}

uint64_t sge2_stream_reader::load_position_file(const char* filename, int32_t icol_strid, int32_t icol_px, int32_t icol_py, bool skip_header, char sep) {
    tsv_reader pos_tr;
    pos_tr.delimiter = (int32_t)sep;
    pos_tr.open(filename);

    if (skip_header) {
        if (!pos_tr.read_line()) {
            error("Cannot read the first line from %s", filename);
        }
        if (pos_tr.nfields <= icol_strid || pos_tr.nfields <= icol_px || pos_tr.nfields <= icol_py) {
            error("The number of fields is : %d, but indices are %d, %d, %d", pos_tr.nfields, icol_strid, icol_px, icol_py);
        }
    }

    return load_position_file_helper(pos_tr, icol_strid, icol_px, icol_py, sep);
}

uint64_t sge2_stream_reader::load_position_file_helper(tsv_reader& pos_tr, int32_t icol_strid, int32_t icol_px, int32_t icol_py, char sep) {
    // need to populate bcd_pos_map
    bcd_pos_map.clear();
    while (pos_tr.read_line()) {
        bcd_pos_map[pos_tr.str_field_at(icol_strid)] = std::make_pair(pos_tr.double_field_at(icol_px), pos_tr.double_field_at(icol_py));
    }
    pos_tr.close();
    return (uint64_t)bcd_pos_map.size();
}

/////////////////////////////////////////////////////////////
// METHODS FOR sge2_stream_writer
/////////////////////////////////////////////////////////////
void sge2_stream_writer::open(const char *bcdf, const char *ftrf, const char *mtxf)
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
    cur_sbcd.nid = nfields = nlines = 0;
}

// close the files
void sge2_stream_writer::close()
{
    if (wh_tmp != NULL)
    {
        if (!flush_mtx())
        {
            error("Cannot flush the mtx file");
        }
    }
    if (wh_ftr != NULL)
    {
        if (hts_close(wh_ftr) != 0)
        {
            error("Cannot close the feature file");
        }
        wh_ftr = NULL;
    }
}

bool sge2_stream_writer::flush_cur_sbcd()
{
    std::string strcnt;

    cat_join_uint64(strcnt, cur_sbcd.cnts, ",");
    if ( (double)cur_sbcd.px == floor((double)cur_sbcd.px) && (double)cur_sbcd.py == floor((double)cur_sbcd.py) )
        hprintf(wh_bcd, "%s\t%.0f\t%.0f\t%s\n", cur_sbcd.strid.c_str(), cur_sbcd.px, cur_sbcd.py, strcnt.c_str());
    else
        hprintf(wh_bcd, "%s\t%.*f\t%.*f\t%s\n", cur_sbcd.strid.c_str(), bcd_precision, cur_sbcd.px, bcd_precision, cur_sbcd.py, strcnt.c_str());

    return true;
}

// add a spatial barcode
bool sge2_stream_writer::add_sbcd(const char *strid, double px, double py) 
{
    if (cur_sbcd.nid > 0)
    {
        flush_cur_sbcd();
    }
    cur_sbcd.strid.assign(strid);
    cur_sbcd.nid++;
    cur_sbcd.px = px;
    cur_sbcd.py = py;

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

bool sge2_stream_writer::write_ftr(const char *id, const char *name, uint64_t nid, std::vector<uint64_t> &cnts)
{
    if (cnts.empty()) {
        cnts.resize(nfields, 0);
    }
    std::string strcnt;
    cat_join_uint64(strcnt, cnts, ",");
    hprintf(wh_ftr, "%s\t%s\t%llu\t%s\n", id, name, nid, strcnt.c_str());
    return true;
}

bool sge2_stream_writer::add_mtx(uint64_t iftr, std::vector<uint64_t> &cnts)
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

bool sge2_stream_writer::add_mtx(uint64_t iftr, std::vector<uint64_t> &cnts, std::vector<int32_t>& icols)
{
    if (nfields == 0)
    {
        nfields = (int32_t)icols.size();
        cur_sbcd.cnts.resize(nfields, 0);
    }
    else if (nfields != (int32_t)icols.size())
    {
        error("The number of fields in the mtx file is not consistent: %d vs %d", nfields, (int32_t)icols.size());
    }

    hprintf(wh_tmp, "%llu %llu", iftr, cur_sbcd.nid);
    for(int32_t i = 0; i < nfields; ++i) {
        if ( icols[i] < 0 || icols[i] >= (int32_t)cnts.size() ) {
            error("Invalid column index %d", icols[i]);
        }
        hprintf(wh_tmp, " %llu", cnts[icols[i]]);
    }
    hprintf(wh_tmp, "\n");

    // std::string strcnt;
    // cat_join_uint64(strcnt, cnts, " ");
    // hprintf(wh_tmp, "%llu %llu %s\n", iftr, cur_sbcd.nid, strcnt.c_str());

    // update ftr_cnts
    if (iftr > ftr_cnts.size() + 1)
        ftr_cnts.resize(iftr);
    if (ftr_cnts[iftr - 1].empty())
        ftr_cnts[iftr - 1].resize(nfields, 0);
    if ((ftr_cnts[iftr - 1].size() != nfields) || (icols.size() != nfields))
        error("The number of fields in the mtx file is not consistent: %d vs %d vs %d", nfields, (int32_t)icols.size(), (int32_t)ftr_cnts[iftr - 1].size());
    for (int32_t i = 0; i < nfields; ++i)
        ftr_cnts[iftr - 1][i] += cnts[icols[i]];

    // update sbcd_cnts
    // assert(cur_sbcd.cnts.size() == nfields);
    for (int32_t i = 0; i < nfields; ++i)
        cur_sbcd.cnts[i] += cnts[icols[i]];

    ++nlines;
    return true;
}

// combine header and contents for mtx to a single compressed mtx.gz file
bool sge2_stream_writer::flush_mtx()
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