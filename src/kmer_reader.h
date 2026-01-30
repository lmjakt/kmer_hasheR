#ifndef _KMER_READER_H
#define _KMER_READER_H

#include <zlib.h>
#include "kseq.h"
#include "suffix_hash.h"

KSEQ_INIT(gzFile, gzread)


// min_ql; a minimum kmer quality expressed as PHRED +33 encoding
// hence, 20, would mean 99% chance of sequence being correct.
/* void kmer_iterator_init(kmer_iterator *it, uint32_t k, unsigned char min_ql); */
/* // these return 1 if a valid k-mer has been found; 0 if past end of sequence */
/* inline int kmer_iterator_begin(kmer_iterator *it, const unsigned char *seq, const unsigned char *qual, uint64_t *kmer_f, uint64_t *kmer_r); */
/* inline int kmer_iterator_nq_begin(kmer_iterator *it, const unsigned char *seq, uint64_t *kmer_f, uint64_t *kmer_r); */
/* inline int kmer_iterator_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r); */
/* inline int kmer_iterator_nq_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r); */

/* typedef enum { SH_single, SH_n } SH_type; */

/* typedef union { */
/*   suffix_hash *sh; */
/*   suffix_hash_n *sh_n; */
/* } SH_union; */


// the pthreads run function will open the file pointed to by file_name
// and use kseq to read sequences
// a kmer will be added to the suffix hash if
// restricted == false || kmer % thread_n == thread_i
// multiple threads will read the same file. This is not elegant, but it seems
// that not much time is spent on the reading and obtaining the k-mers, so this
// may be still scale OK. We'll see.
typedef struct {
  const char *file_name;
  //  gzFile fp;
  //  kseq_t *kseq;
  int k;
  uint32_t suffix_bits;
  unsigned char min_q;
  unsigned char restricted;
  uint32_t thread_n;
  uint32_t thread_i;
  pthread_t thread;
  uint32_t source;
  suffix_hash_n *sh;
  size_t max_reads;
} kmer_reader;

void init_kmer_reader(kmer_reader *kr, const char *file, suffix_hash_n *sh,
		      int k, uint32_t suffix_bits, uint32_t source,
		      unsigned char min_q, unsigned char restricted, uint32_t thread_n,
		      uint32_t thread_i, size_t max_reads);

void *kmer_reader_read(void *v_args);

// a collection of readers
typedef struct {
  uint32_t thread_n;
  kmer_reader *readers;
  // if thread_n = 1, then restricted -> false
} kmer_reader_pool;

// this will start the threads running; we need a function to join them as well
suffix_hash_n *init_kmer_reader_pool(kmer_reader_pool *krp, const char *file_name,
				     int k, uint32_t prefix_bits, size_t max_size, uint32_t thread_n,
				     unsigned char min_q, size_t max_reads,
				     uint32_t source_n, uint32_t source);

// This takes an existing suffix hash and extends it.
suffix_hash_n* init_kmer_reader_pool_sh(kmer_reader_pool *krp, const char *file_name,
					int k, suffix_hash_n *sh, size_t max_size, uint32_t thread_n,
					unsigned char min_q, size_t max_reads, uint32_t source);


void join_kmer_reader_pool(kmer_reader_pool *krp);
void free_kmer_reader_pool(kmer_reader_pool *krp);


// return the counts of the kmers at all positions; 
int seq_kmer_counts(const char* seq, size_t seq_l, int* counts, suffix_hash_n *sh, int k);


#endif
