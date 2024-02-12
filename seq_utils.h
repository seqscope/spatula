#ifndef __SEQ_UTILS_H
#define __SEQ_UTILS_H

#include <cstdint>

int32_t seq_iupac_mismatch(const char *seq, const char *pattern, int32_t len);
void seq_revcomp(char *seq, int32_t len);
uint64_t seq2bits(const char *seq, int32_t len, uint8_t nonACGTs = 0x03);
uint64_t seq2nt5(const char *seq, int32_t len);

#endif // __SEQ_UTILS_H
