#include "Q_to_log_likelihood.h"
#include "kmer_util.h"

size_t skip_n(const char *seq, size_t i){
  while(seq[i] && LC(seq[i]) == 'n')
    ++i;
  return(i);
}

size_t skip_n_qual(const char *seq, const char *qual, char min_q, size_t i){
  while(seq[i] && (LC(seq[i]) == 'n' || (qual && qual[i] < min_q)))
    ++i;
  return(i);
}

// returns the value of j; if it runs out of sequence (i.e. last k-mer)
// or it finds an N. The caller must check the return value.
size_t init_kmer(const char *seq, size_t i, unsigned long *offset, int k){
  size_t j = 0;
  while(seq[i]){
    *offset = 0;
    for(j=0; j < k && seq[i+j] && LC(seq[i+j]) != 'n'; ++j){
      *offset = UPDATE_OFFSET(*offset, seq[i+j]);
    }
    if(seq[i+j] == 0 || j == k)
      break;
    // otherwise we hit Ns again;
    i = skip_n(seq, i + j);
    j=0;
  }
  return( i + j );
}

// do forward and reverse complement
size_t init_kmer_qual_2(const char *seq, const char *qual, char min_q, size_t i,
			uint64_t *offset, uint64_t *offset_rc,
			int k){
  size_t j = 0;
  while(seq[i]){
    *offset = 0;
    *offset_rc = 0;
    for(j=0; j < k && seq[i+j] && LC(seq[i+j]) != 'n' && (!qual || qual[i+j] >= min_q); ++j){
      *offset = UPDATE_OFFSET(*offset, seq[i+j]);
      *offset_rc = UPDATE_OFFSET_RC(*offset_rc, seq[i+j]);
    }
    if(seq[i+j] == 0 || j == k) // 
      break;
    // otherwise we hit Ns again;
    i = skip_n_qual(seq, qual, min_q, i + j);
    j = 0;
  }
  return( i + j );
}

void kmer_iterator_init(kmer_iterator *it, uint32_t k, unsigned char min_ql){
  it->seq = 0; it->qual = 0;
  it->fwd = 0; it->rev = 0;
  it->k = k;
  it->kmer_ll = 0;
  it->prev_ll = 0;
  it->min_ll = q_to_ll[min_ql];
  it->kmer_mask = (1ULL << (k*2)) - 1;
  it->rc_shift = 64 - (k * 2);
}

// this is called from kmer_interator_begin when qual == 0
// or from kmer_iterator_nq_next;
int kmer_iterator_nq_begin(kmer_iterator *it, const unsigned char* seq, uint64_t *kmer_f, uint64_t *kmer_r){
  it->fwd = 0;  it->rev = 0;
  it->seq = seq;
  size_t i = 0;
  while(*seq && LC(*seq) != 'n' && i < it->k){
    it->fwd = UPDATE_OFFSET(it->fwd, *seq);
    it->rev = UPDATE_OFFSET_RC(it->rev, *seq);
    ++seq; ++i;
  }
  // if i == k, then set the it->seq and it->qual
  if(i == it->k){
    it->seq = seq;
    *kmer_f = it->fwd & it->kmer_mask;
    *kmer_r = (it->rev >> it->rc_shift);
    return(1);
  }
  // else, skip low quality bases and recurse:
  while(*seq && LC(*seq) == 'n')
    ++seq;

  if(!*seq)
    return(0);
  // possibly consider writing this as a loop to avoid deep recursion and
  // stack overflow; 
  return(kmer_iterator_nq_begin(it, seq, kmer_f, kmer_r));
}

int kmer_iterator_begin(kmer_iterator *it, const unsigned char* seq, const unsigned char* qual, uint64_t *kmer_f, uint64_t *kmer_r){
  if(!qual)
    return(kmer_iterator_nq_begin(it, seq, kmer_f, kmer_r));
  
  it->fwd = 0;  it->rev = 0;
  it->seq = seq; it->qual = qual;
  size_t i = 0;
  double kmer_ll = 0;
  double prev_ll = 0;
  while(*seq && ((kmer_ll = (kmer_ll + q_to_ll[*qual])) > it->min_ll) && i < it->k){
    it->fwd = UPDATE_OFFSET(it->fwd, *seq);
    it->rev = UPDATE_OFFSET_RC(it->rev, *seq);
    prev_ll = q_to_ll[*qual];
    ++seq; ++qual; ++i;
  }
  // if i == k, then set the it->seq and it->qual
  if(i == it->k){
    it->seq = seq; it->qual = qual;
    it->prev_ll = prev_ll;
    it->kmer_ll = kmer_ll;
    *kmer_f = it->fwd & it->kmer_mask;
    *kmer_r = (it->rev >> it->rc_shift);
    return(1);
  }
  // else, skip low quality bases and recurse:
  while(*seq && q_to_ll[*qual] <= it->min_ll){
    ++seq; ++qual;
  }
  if(!*seq)
    return(0);
  // possibly consider writing this as a loop to avoid deep recursion and
  // stack overflow; 
  return(kmer_iterator_begin(it, seq, qual, kmer_f, kmer_r));
}


int kmer_iterator_nq_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r){
  if(!*it->seq)
    return(0);
  if(LC(*it->seq) == 'n')
    return(kmer_iterator_nq_begin(it, ++(it->seq), kmer_f, kmer_r));
  // update the iterator and return 1;
  it->fwd = UPDATE_OFFSET(it->fwd, *it->seq);
  it->rev = UPDATE_OFFSET_RC(it->rev, *it->seq);
  *kmer_f = it->fwd & it->kmer_mask;
  *kmer_r = it->rev >> it->rc_shift;
  ++(it->seq);
  return(1);
}

int kmer_iterator_next(kmer_iterator *it, uint64_t *kmer_f, uint64_t *kmer_r){
  if(!it->qual)
    return(kmer_iterator_nq_next(it, kmer_f, kmer_r));
  if(!*it->seq)
    return(0);
  it->kmer_ll += (q_to_ll[*it->qual] - it->prev_ll);
  if(it->kmer_ll < it->min_ll)
    return(kmer_iterator_begin(it, ++(it->seq), ++(it->qual), kmer_f, kmer_r));
  // update the iterator and return 1;
  it->fwd = UPDATE_OFFSET(it->fwd, *it->seq);
  it->rev = UPDATE_OFFSET_RC(it->rev, *it->seq);
  it->prev_ll = q_to_ll[*it->qual];
  *kmer_f = it->fwd & it->kmer_mask;
  *kmer_r = it->rev >> it->rc_shift;
  it->seq++; it->qual++;
  return(1);
}

