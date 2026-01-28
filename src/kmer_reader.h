#ifndef _KMER_READER_H
#define _KMER_READER_H

#include <zlib.h>
#include "kseq.h"
#include "suffix_hash.h"

KSEQ_INIT(gzFile, gzread)

// a macro to update an offset as defined here
#define UPDATE_OFFSET(off, c) ( ((off) << 2) | (( (c) >> 1) & 3) )
#define UPDATE_OFFSET_RC(off, c) ( ((off) >> 2) | ( (( (((uint64_t)(c) >> 1) & 3) + 2) % 4) << 62 ))
#define LC(c) ( (c) | 0x20 )
#define UC(c) ( (c) & 0xCF )
#define MAX_K 32

  
size_t skip_n(const char *seq, size_t i);
size_t skip_n_qual(const char *seq, const char *qual, char min_q, size_t i);

size_t init_kmer(const char *seq, size_t i, unsigned long *offset, int k);

// do forward and reverse complement
size_t init_kmer_qual_2(const char *seq, const char *qual, char min_q, size_t i,
			uint64_t *offset, uint64_t *offset_rc,
			int k);

// to iterate over a sequence emitting k-mers we define an iterator struct
// seq and qual: 0 terminated strings giving sequence and qualities
//               both must be defined and be of equal length.
// i           : the current position in the sequence (the last position of the base)
// fwd         : the current full sized k-mer (unmasked for consistency with rev)
// rev         : the kmer for the reverse complement: this needs to be shifted and masked
//               to give the correct k-mer
// k           : the kmer -length
// kmer_b_ll  : for each base in the sequence, the log likelihood of a base being correct
//               kmer_b_ll[i] = log(1-10^(-Q[i]/10)). The sum of this will give the estimated
//               likelihood that the kmer is correct, and used to skip bases.
//               This will be used as a circular buffer with values added at (i % k)
// kmer_ll : the sum of kmer_b_ll, i.e. the log likelihood this kmer sequence is correct
//                this will need to be updated for each new base.
//                update by subtracting the old value at kmer_b_ll(i%k) and adding the new value
// q_to_cor_lp  : NOTE, defined statically in the c file.
//                an array for converting phread quality scores (phred33 encoding assumed)
//                to log likelihoods of being correct; for '!' and lower this will be set to
//                log(DBL_MIN). It might be speed things up if this is defined statically
//                but for now, simply define when instantiated.
// min_cor_ll    : the minimum log likelihood of being correct for including a k-mer in the count

typedef struct {
  const unsigned char *seq;
  const unsigned char *qual;
  uint64_t fwd;
  uint64_t rev;
  uint32_t k;
  double prev_ll;
  double kmer_ll;
  double min_ll;
  uint64_t kmer_mask;
  uint32_t rc_shift;
  // q_to_cor_p[256]
  // is defined globally in the c file.
} kmer_iterator;

// min_ql; a minimum kmer quality expressed as PHRED +33 encoding
// hence, 20, would mean 99% chance of sequence being correct.
/* void kmer_iterator_init(kmer_iterator *it, uint32_t k, unsigned char min_ql); */
/* // these return 1 if a valid k-mer has been found; 0 if past end of sequence */
/* inline int kmer_iterator_begin(kmer_iterator *it, const unsigned char *seq, const unsigned char *qual, uint64_t *kmer_f, uint64_t *kmer_r); */
/* inline int kmer_iterator_nq_begin(kmer_iterator *it, const unsigned char *seq, uint64_t *kmer_f, uint64_t *kmer_r); */
/* inline int kmer_iterator_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r); */
/* inline int kmer_iterator_nq_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r); */

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


// return the counts of the kmers at all positions; 
int seq_kmer_counts(const char* seq, size_t seq_l, int* counts, suffix_hash *sh, int k);


#endif
