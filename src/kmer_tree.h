#ifndef _KMER_TREE_H
#define _KMER_TREE_H

#include <stdint.h>

// This defines a two level hash for holding counts of kmers
// The data structure is similar to the one used by the program
// meryl.
// It allows direct access to counts through two memory address
// hops.

#define MAX_PRE_BITS 36
#define MAX_SUF_BITS 36
#define KMAX_BITS = 64;


typedef struct {
  uint32_t suffix_bits;
  size_t suffix_n;
  uint32_t *suffix_counts;
} suf_counts;

suf_counts *init_suf_counts(uint32_t bits);
void free_suf_counts(suf_counts *counts);

// a should be a kmer_tree data structure
// this gives the amount of memory currently used
#define KTREE_MEM(a) ( 2^((a).prefix_bits) * 8 + 2^((a).suffix_bits) * 4 * (a).allocated )
// This gives the amount of memory required if one more suffix array is allocated
#define KTREE_MEM_A(a) ( 2^((a).prefix_bits) * 8 + 2^((a).suffix_bits) * 4 * (1 + (a).allocated) )

typedef struct {
  uint32_t suffix_bits;
  int max_count;
  uint64_t max_count_kmer;
  uint64_t kmer_mask;
  uint64_t suffix_mask;
  uint32_t prefix_bits;
  suf_counts **prefixes;
  size_t prefix_n;
  size_t allocated;
  size_t max_size;
} kmer_tree;

kmer_tree init_kmer_tree(uint32_t prefix_bits, uint32_t suffix_bits, size_t max_size);
void free_kmer_tree(kmer_tree *kt);
int add_kmer(kmer_tree *kt, uint64_t kmer);

// counts should point to a an array of integers of length counts_n
// counts will not be cleared by this function
int kmer_count(kmer_tree *kt, uint64_t kmer);
size_t count_spectrum(kmer_tree *kt, double *counts, uint32_t counts_n);

#endif
