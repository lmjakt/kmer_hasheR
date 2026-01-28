#ifndef _SUFFIX_HASH_H
#define _SUFFIX_HASH_H

#include <stdint.h>
#include "khash.h"
#include "thread_queue.h"

// This defines a two level hash for holding counts of kmers
// The data structure is similar to the one used by the program
// meryl.
// It allows direct access to counts through two memory address
// hops.

#define SH_MAX_PRE_BITS 36
#define SH_MAX_SUF_BITS 32
#define SH_KMAX_BITS 64

// like kcount, but hold two two uint32_t counts
// This is ugly as hell, but it should reduce the
// amount of pointers. 
typedef struct {
  uint32_t n[2];
} count_2_t;

typedef struct {
  uint32_t n[3];
} count_3_t;

typedef struct {
  uint32_t n[4];
} count_4_t;


// each prefix should point to an int int hash:
KHASH_MAP_INIT_INT(kcount, uint32_t)
KHASH_MAP_INIT_INT(kcount_2, count_2_t)
KHASH_MAP_INIT_INT(kcount_3, count_3_t)
KHASH_MAP_INIT_INT(kcount_4, count_4_t)

typedef struct {
  uint32_t suffix_bits;
  int max_count;
  uint64_t max_count_kmer;
  uint64_t kmer_mask;
  uint64_t suffix_mask;
  uint32_t prefix_bits;
  khash_t(kcount) **prefixes;
  size_t prefix_n;
  size_t allocated;
  size_t max_size;
} suffix_hash;

typedef struct {
  uint32_t k;
  uint32_t prefix_bits;
  uint32_t suffix_bits;
  uint32_t counts_n; // can be, 1, 2, 3, or 4
  uint64_t kmer_mask;
  uint64_t suffix_mask;
  size_t prefix_n;
  void **prefixes;  // khash_t(kcount[2|3|4])*
  uint64_t *kmer_counts;
} suffix_hash_n;
  
// use for merging operations like union intersect set diff and so on.
/* typedef struct { */
/*   suffix_hash *sh1; */
/*   suffix_hash *sh2; */
/*   uint32_t thread_i; */
/*   uint32_t thread_n; */
/*   pthread_t thread; */
/*   khash_t(kcount_2) **prefixes; // write counts here..  */
/* } sh_pair_t; */

suffix_hash init_suffix_hash(uint32_t prefix_bits, uint32_t suffix_bits, size_t max_size);
void init_suffix_hash_p(suffix_hash *sh, uint32_t prefix_bits, uint32_t suffix_bits, size_t max_size);
void free_suffix_hash(suffix_hash *kt);
int sh_add_kmer(suffix_hash *sh, uint64_t kmer);

uint32_t sh_kmer_count(suffix_hash *sh, uint64_t kmer);
// counts should point to a an array of integers of length counts_n
size_t sh_count_spectrum(suffix_hash *sh, double *counts, uint32_t counts_n);


int init_suffix_hash_n(suffix_hash_n* sh, uint32_t counts_n,
			uint32_t k, uint32_t prefix_bits, uint32_t suffix_bits);
void free_suffix_hash_n(suffix_hash_n *sh);
int sh_n_add_kmer(suffix_hash_n *sh, uint32_t source, uint64_t kmer);

// counts must point to an array of length sh->counts_n
void sh_kmer_count_n(suffix_hash_n *sh, uint64_t kmer, uint32_t* counts);
// counts should be an array of length counts_l * sh->counts_n
// comb defines which k-mers are counted. For a kmer to be counted, it must have
// been observed in all bits set to 1, and NOT observed in any bit set to 0
// if inner is true
// if inner is false, then all k-mers counted in any of the bits will be included.
size_t sh_count_spectrum_nc(suffix_hash_n *sh, uint32_t *counts, uint32_t counts_l,
			    uint32_t comb_in, uint32_t inner);



  
/* typedef struct { */
/*   sh_pair_t *threads; */
/*   uint32_t thread_n; */
/* } sh_set_op_t; */

//sh_set_op_t init_sh_set_op_t(suffix_hash_2* shp, uint32_t thread_n, typeof(void *(void *)) *start_routine);




// Set operations on pairs of suffix_hashes
// args: a pointer to a sh_pair_t
//void *sh_intersect_thread( void* args );
//suffix_hash_2 sh_intersect(suffix_hash *sh1, suffix_hash *sh2, uint32_t thread_n);
//

// structs and functions for using pthreads
// The arguments that should be given to threads that will count
// kmers.
// WARNING: this design assumes that each thread has a distinct set
// of hashes that it owns. But it is up to the reading thread to pass
// the appropriate k-mers to this structure.
typedef struct {
  uint32_t prefix_bits;
  uint32_t suffix_bits;
  uint64_t suffix_mask;
  khash_t(kcount) **hashes;
  kmer_queue *queue;
  size_t hashes_allocated;
  size_t kmer_n;
  int *done; // if not 0 exit thread.
  size_t thread_i;
} suffix_hash_t_args;

void init_thread_args(suffix_hash_t_args *args,
		      uint32_t prefix_bits, uint32_t suffix_bits,
		      khash_t(kcount) **hashes, kmer_queue *queue,
		      int *done, size_t thread_i);

// add k-mers in the queue (i.e. add to the suffix array) in a separate thread
// call pthread_create( ) with this as an argument;
// args should be a suffix_hash_t_args structure
void *sh_process_kmers(void *args);


typedef struct {
  uint32_t suffix_bits;
  uint32_t prefix_bits;
  uint64_t kmer_mask;
  size_t prefix_n; // the number of hashes
  khash_t(kcount) **hashes;
  size_t nthreads;
  suffix_hash_t_args *t_args;
  pthread_t *threads;
  kmer_queue *queues;
  int done;
} suffix_hash_mt;

void init_suffix_hash_mt(suffix_hash_mt* sh,
			 uint32_t prefix_bits, uint32_t suffix_bits,
			 size_t nthreads, size_t queue_buffer_size);

void free_suffix_hash_mt(suffix_hash_mt *sh);

int suffix_hash_mt_add(suffix_hash_mt* sh, uint64_t);
// void suffix_hash_mt_join(suffix_hash_mt* sh);

#endif
