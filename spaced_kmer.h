#ifndef __SPACED_KMER_H
#define __SPACED_KMER_H

#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>
#include <functional>
#include <stdint.h>

#include "Error.h"

class array_of_bytes_t {
public:
  int32_t nbytes;     // byte size for individual units
  int32_t chunk_size; // bytes per chunk
  int64_t nelems;     // current number of elements
  std::vector<void*> chunks; // chunks to store elements

  array_of_bytes_t(int32_t _nbytes, int32_t _chunk_size) : nbytes(_nbytes), chunk_size(_chunk_size), nelems(0) {
    if ( chunk_size % nbytes != 0 )
      error("[E:%s:%d %s] chunk_size = %d is not divisible by nbytes = %d", __LINE__, __FILE__, __PRETTY_FUNCTION__, chunk_size, nbytes);
  }
  
  ~array_of_bytes_t() {
    for(std::vector<void*>::iterator it = chunks.begin(); it != chunks.end(); ++it) {
      free(*it);
    }
  }

  inline uint64_t add_elem(void* elem) {
    if ( ( nelems * nbytes ) % chunk_size == 0 ) { // add a new chunk
      void* p = malloc(chunk_size);
      chunks.push_back(p);
    }
    void* dest = (void*)((char*)chunks.back() + (nelems * nbytes) % chunk_size);
    memcpy(dest, elem, nbytes);
    return nelems++;
  }

  inline void* get_elem(uint64_t idx) {
    return (void*)((char*)chunks[(idx * nbytes) / chunk_size] + ( idx * nbytes ) % chunk_size);
  }
};

typedef std::unordered_map<uint64_t, std::vector<uint64_t> > hashtable_t;
typedef std::unordered_map<uint64_t, std::vector<uint64_t> >::iterator hashtable_it_t;

class spaced_kmer_t {
public:
  int32_t lseq;   // length of sequence (e.g. 32)
  int32_t lunit;  // length of kmer as a 'unit' of sequence (e.g. 8)
  int32_t nunit;  // number of units in each kmer entry (e.g. 2)
  std::vector< std::vector<int32_t> > conf; // configurations of units
  std::vector<hashtable_t> tables; // tables of k-mers for each configuration
  array_of_bytes_t seq_2bits; // 2 bit sequences, 4 nucleotides per byte

  // initialize configuration
  spaced_kmer_t(int32_t _lseq, int32_t _lunit, int32_t _nunit, int32_t _chunk_size = 100000);

  uint64_t masked_key(void* p2bit, const std::vector<int32_t>& mask);
  uint64_t add_seq(const char* seq);
  int32_t query_seq(const char* seq, std::set<uint64_t>& rst);

  const char* get_seq_by_index(uint64_t idx, char* dst = NULL);
  int32_t conv_seq_to_2bit(const char* seq, void* dst = NULL);
  int32_t conv_2bit_to_seq(const void* p2bit, char* dst = NULL);

private:
  char buf[255];
  int32_t nuc2bit[255];
  char acgts[4];
  std::hash<std::string> str_hash;  
};

#endif // SPACED_KMER_H
