#ifndef __FASTQ_UTILS_H
#define __FASTQ_UTILS_H

#include <cstdint>
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"

struct _fq_rec
{
    kstring_t rname;
    kstring_t seq;
    kstring_t comment;
    kstring_t qual;

    // default constructor
    _fq_rec()
    {
        rname.l = 0;
        rname.m = 0;
        rname.s = NULL;

        seq.l = 0;
        seq.m = 0;
        seq.s = NULL;

        comment.l = 0;
        comment.m = 0;
        comment.s = NULL;

        qual.l = 0;
        qual.m = 0;
        qual.s = NULL;
    }

    // destructor, free the memory
    ~_fq_rec()
    {
        if (rname.s != NULL)
            free(rname.s);
        if (seq.s != NULL)
            free(seq.s);
        if (comment.s != NULL)
            free(comment.s);
        if (qual.s != NULL)
            free(qual.s);
    }
};
typedef struct _fq_rec fq_rec_t;

bool read_fastq_record(htsFile *hf, fq_rec_t &rec);
bool write_fastq_record(htsFile *hf, fq_rec_t &rec);

#endif // __FASTQ_UTILS_H
