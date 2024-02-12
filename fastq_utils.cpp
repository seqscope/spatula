#include "seq_utils.h"
#include "fastq_utils.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include <cstdio>
#include <cstring>

// Read a record from a FASTQ file
bool read_fastq_record(htsFile *hf, fq_rec_t &rec)
{
//    notice("foo");
    int32_t len = hts_getline(hf, KS_SEP_LINE, &rec.rname);
    if (len <= 0)
    { // EOF reached. No record is returned
        return false;
    }
//    notice("bar %d", len);
    len = hts_getline(hf, KS_SEP_LINE, &rec.seq);
    if (len == 0)
        error("Unexpected EOF in reading a FASTQ file");
    len = hts_getline(hf, KS_SEP_LINE, &rec.comment);
    if (len == 0)
        error("Unexpected EOF in reading a FASTQ file");
    if (rec.comment.s[0] != '+')
        error("Unexpected line in reading a FASTQ file. Expected '+', but observed %s.", rec.qual.s);
    len = hts_getline(hf, KS_SEP_LINE, &rec.qual);
    if (len == 0)
        error("Unexpected EOF in reading a FASTQ file");
//    notice("baz %d", len);
    return true;
}

bool write_fastq_record(htsFile* wf, fq_rec_t & rec)
{
    if ( wf != NULL ) {
        hprintf(wf, "%s\n", rec.rname.s);
        hprintf(wf, "%s\n", rec.seq.s);
        hprintf(wf, "%s\n", rec.comment.s);
        hprintf(wf, "%s\n", rec.qual.s);
        return true;
    }
    else return false;
}