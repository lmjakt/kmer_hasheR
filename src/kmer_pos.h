#ifndef _KMER_POS_H
#define _KMER_POS_H

#include <stdint.h>
#include "khash.h"
#include "kvec.h"

#define OPT_KMER 0
#define OPT_POS 1
#define OPT_PAIRS 2
#define OPT_COUNT 3
#define N_OPTS 4

// The following are not currently used in kmer_pos.c.
// Commented and moved to kmer_hash.c to avoid compiler
// warnings.
// Consider moving functions to kmer_pos.c in the future
/* static const uint32_t pos_opt_flags[N_OPTS] = {1, 2, 4, 8}; */
/* static const char* kmer_pos_fields[N_OPTS] = {"kmer", "pos", "pair.pos", "count"}; */


// This will translate into code that defines an
// hash using 64 bit integers as keys
// and the value as a kmer_pos_t
// It is intended for relatively short sequences.
// kvec_t(int)
// the flag value can be set by functions in order to
// mark k-mers that should be treated differently.
typedef struct {
  kvec_t(int) v;
  uint64_t kmer;
  uint64_t kmer_flag;
} kmer_pos_t;

KHASH_MAP_INIT_INT64(kmer_h, kmer_pos_t)


typedef struct {
  khash_t(kmer_h) *hash;
  int k;
  size_t kmer_count;
  int sorted;
} khash_ptr;

void clear_kmer_h(khash_t(kmer_h) *hash);

void sort_kmer_pos(khash_ptr *hash_ptr);

int seq_to_hash(const char *seq, int k, khash_t(kmer_h) *hash);

#endif
