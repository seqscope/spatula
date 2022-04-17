#include <algorithm>
#include <cassert>
#include <cstring>
#include <set>
#include "spaced_kmer.h"

spaced_kmer_t::spaced_kmer_t(int32_t _lseq, int32_t _lunit, int32_t _nunit, int32_t _chunk_size) : lseq(_lseq), lunit(_lunit), nunit(_nunit), seq_2bits((_lseq+3)/4, _chunk_size) {
  // lseq - total sequence length            (e.g. 32, 20)
  // lunit - smallest unit for spaced-kmer   (e.g. 8,  4)
  // nunit - number of smallest unit to make up a spaced kmer (e.g. 2, 3)
  // number of configuration == (lseq/lunit) choose (nunit)
  int32_t n_total_units = (lseq + (lunit-1))/ lunit;

  if ( lunit % 4 != 0 )
    error("[E:%s:%d %s], Invalid lunit=%d. Currently lunits multiples of 4 are supported only", __FILE__, __LINE__, __PRETTY_FUNCTION__, lunit);

  if ( lunit * nunit > 32 )
    error("[E:%s:%d %s], Invalid (lunit,nunit)=(%d,%d). Currently lunits*nunit cannot be greater than 32", __FILE__, __LINE__, __PRETTY_FUNCTION__, lunit, nunit);
 

  std::string bitmask(nunit, 1); // K leading 1's
  bitmask.resize(n_total_units, 0); // N-K trailing 0's
  
  // print integers and permute bitmask
  do {
    std::vector<int32_t> comb;
    for (int32_t i = 0; i < n_total_units; ++i) { // [0..N-1] integers
	if (bitmask[i]) comb.push_back(1);
	else comb.push_back(0);
    }
    conf.push_back(comb);
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  tables.resize(conf.size());

  // initialize nucbit
  for(int32_t i=0; i < 255; ++i) nuc2bit[i] = -1;
  nuc2bit['A'] = 0; nuc2bit['a'] = 0;
  nuc2bit['C'] = 1; nuc2bit['c'] = 1;
  nuc2bit['G'] = 2; nuc2bit['g'] = 2;
  nuc2bit['T'] = 3; nuc2bit['t'] = 3;

  acgts[0] = 'A'; acgts[1] = 'C'; acgts[2] = 'G'; acgts[3] = 'T';

  //notice("nunit = %d, n_total_units = %d", nunit, n_total_units);
  //notice("conf.size() = %u", conf.size());
  //notice("conf[0].size() = %u", conf[0].size());  
  //notice("conf[0] = (%d,%d,%d,%d,%d)", conf[0][0], conf[0][1], conf[0][2], conf[0][3],conf[0][4]);  
}

int32_t spaced_kmer_t::conv_2bit_to_seq(const void* p2bit, char* dst) {
  const char* psrc = (const char*)p2bit;
  char* pdst = dst != NULL ? dst : buf;
  int32_t i;
  for(i=0; i < lseq; ++i) {
    pdst[i] = acgts[ (psrc[i/4] >> (6 - (i%4)*2)) & 0x03 ];
  }
  pdst[lseq] = '\0';
  return i;
}

int32_t spaced_kmer_t::conv_seq_to_2bit(const char* seq, void* dst) {
  const char* s = seq;
  char* pdst = dst != NULL ? (char*)dst : buf;
  int32_t i, twobit;
  for(i=0; *s != '\0'; ++i, ++s) {
    twobit = nuc2bit[*s];
    if ( twobit < 0 )
      twobit = str_hash(s) % 4;
    if ( i % 4 == 0 ) pdst[i/4] = (twobit << 6);
    else pdst[i/4] |= (twobit << (6 - (i%4)*2 ));
  }
  return i;
}

uint64_t spaced_kmer_t::masked_key(void* p2bit, const std::vector<int32_t>& mask) {
  int32_t bytePerUnit = lunit/4;
  uint64_t key = 0;
  const char* s = (const char*)p2bit;
  for(int32_t j=0, k=0; j < (int32_t)mask.size(); ++j) {
    if ( mask[j] == 1 ) {
      for(int32_t b=0; b < bytePerUnit; ++b) {
	key = (key << 8) | (s[bytePerUnit*j + b] & 0x00ff);
      }
      //fprintf(stderr,"%d %d %llx\n", j, k, key);
      ++k;
    }
  }
  return key;
}

// add A/C/G/T sequences to build spaced kmers
uint64_t spaced_kmer_t::add_seq(const char* seq) {
  if ( conv_seq_to_2bit(seq) != lseq )
    error("[E:%s:%d %s] seq=%s does not have expected length %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, seq, lseq);

  uint64_t idx = seq_2bits.add_elem(buf);
  void* p2bit = seq_2bits.get_elem(idx);
  
  int32_t bytePerUnit = lunit/4;
  for(int32_t i=0; i < (int32_t)conf.size(); ++i) {
    uint64_t key = masked_key(p2bit, conf[i]);
    //notice("%s\t%u %u %u %u %u\t%llx\t%llu", seq, (uint8_t)buf[0], (uint8_t)buf[1], (uint8_t)buf[2], (uint8_t)buf[3], (uint8_t)buf[4], key, key);
    tables[i][key].push_back(idx);
  }
  //if ( rand() % 20 == 0 ) error("halt");
}

const char* spaced_kmer_t::get_seq_by_index(uint64_t idx, char* dst) {
  char* pdst = dst != NULL ? dst : buf;
  void* p2bit = seq_2bits.get_elem(idx);
  conv_2bit_to_seq(p2bit, pdst);
  return pdst;
}

int32_t spaced_kmer_t::query_seq(const char* seq, std::set<uint64_t>& rst) {
  char pbuf[255];
  if ( conv_seq_to_2bit(seq, (void*)pbuf) != lseq )
    error("[E:%s:%d %s] seq=%s does not have expected length %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, seq, lseq);

  int32_t bytePerUnit = lunit/4;
  for(int32_t i=0; i < (int32_t)conf.size(); ++i) {
    uint64_t key = masked_key((void*)pbuf, conf[i]);
    hashtable_it_t itr = tables[i].find(key);
    if ( itr != tables[i].end() ) { // found something
      //notice("i = %d, key = %x, found = %zu", i, key, itr->second.size());
      for(std::vector<uint64_t>::iterator it = itr->second.begin();
	  it != itr->second.end(); ++it) {
	rst.insert(*it);
      }
    }
  }
  return (int)rst.size();
}
