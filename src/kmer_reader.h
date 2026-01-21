#ifndef _KMER_READER_H
#define _KMER_READER_H

#include "kseq.h"
#include "suffix_hash.h"

KSEQ_INIT(gzFile, gzread)

// a macro to update an offset as defined here
#define UPDATE_OFFSET(off, c) ( ((off) << 2) | (( (c) >> 1) & 3) )
#define UPDATE_OFFSET_RC(off, c) ( ((off) >> 2) | ( (( (((uint64_t)(c) >> 1) & 3) + 2) % 4) << 62 ))
#define LC(c) ( (c) | 0x20 )
#define UC(c) ( (c) & 0xCF )
#define MAX_K 32

size_t skip_n_qual(const char *seq, const char *qual, char min_q, size_t i);

// do forward and reverse complement
size_t init_kmer_qual_2(const char *seq, const char *qual, char min_q, size_t i,
			uint64_t *offset, uint64_t *offset_rc,
			int k);

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
  suffix_hash *sh;
  size_t max_reads;
} kmer_reader;

void init_kmer_reader(kmer_reader *kr, const char *file, suffix_hash *sh, int k, uint32_t suffix_bits,
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
suffix_hash *init_kmer_reader_pool(kmer_reader_pool *krp, const char *file_name,
				   int k, uint32_t prefix_bits, size_t max_size, uint32_t thread_n,
				   unsigned char min_q, size_t max_reads);

// This takes an existing suffix hash and extends it.
suffix_hash *init_kmer_reader_pool_sh(kmer_reader_pool *krp, const char *file_name,
				      int k, suffix_hash *sh, size_t max_size, uint32_t thread_n,
				      unsigned char min_q, size_t max_reads);


void join_kmer_reader_pool(kmer_reader_pool *krp);
void free_kmer_reader_pool(kmer_reader_pool *krp);
  
#endif
