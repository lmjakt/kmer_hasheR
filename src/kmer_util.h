#ifndef _KMER_MACROS_H
#define _KMER_MACROS_H

#include <stdint.h>
#include <stddef.h>

// macros that update 64 bit encodings of 32-mers
#define UPDATE_OFFSET(off, c) ( ((off) << 2) | (( (c) >> 1) & 3) )
#define UPDATE_OFFSET_RC(off, c) ( ((off) >> 2) | ( (( (((uint64_t)(c) >> 1) & 3) + 2) % 4) << 62 ))
#define LC(c) ( (c) | 0x20 )
#define UC(c) ( (c) & 0xCF )
#define MAX_K 32

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


size_t skip_n(const char *seq, size_t i);
size_t skip_n_qual(const char *seq, const char *qual, char min_q, size_t i);

size_t init_kmer(const char *seq, size_t i, unsigned long *offset, int k);

// do forward and reverse complement
size_t init_kmer_qual_2(const char *seq, const char *qual, char min_q, size_t i,
			uint64_t *offset, uint64_t *offset_rc,
			int k);

// k-mer iterator functions; for iterating a k-mer along a sequence
// rejecting kmer-s with N, or k-mers with too large a probability of
// containing errors:
void kmer_iterator_init(kmer_iterator *it, uint32_t k, unsigned char min_ql);

int kmer_iterator_nq_begin(kmer_iterator *it, const unsigned char* seq,
			   uint64_t *kmer_f, uint64_t *kmer_r);

int kmer_iterator_begin(kmer_iterator *it, const unsigned char* seq,
			const unsigned char* qual, uint64_t *kmer_f, uint64_t *kmer_r);

int kmer_iterator_nq_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r);

int kmer_iterator_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r);


#endif
