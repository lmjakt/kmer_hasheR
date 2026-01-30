#include <stdio.h>
#include <zlib.h>
#include <float.h> // for DBL_MAX
#include <math.h>
#include "kmer_util.h"
#include "kmer_reader.h"

void init_kmer_reader(kmer_reader *kr, const char *file, suffix_hash_n *sh,
		      int k, uint32_t suffix_bits, uint32_t source,
		      unsigned char min_q, unsigned char restricted, uint32_t thread_n,
		      uint32_t thread_i, size_t max_reads){
  kr->file_name = file;
  kr->sh = sh;
  kr->source = source;
  //  kr->sh = sh;
  kr->k = k;
  kr->suffix_bits = suffix_bits;
  kr->min_q = min_q;
  kr->restricted = restricted;
  kr->thread_n = thread_n;
  kr->thread_i = thread_i;
  kr->max_reads = max_reads;
  // the pthread_t thread is initialised when pthread_create is called
}

int kr_add_kmer(kmer_reader *kr, uint64_t kmer_f, uint64_t kmer_r,
		size_t *word_count, size_t *kmer_count, uint32_t source){
  uint64_t kmer = kmer_f < kmer_r ? kmer_f : kmer_r;
  size_t prefix_i = kmer >> kr->suffix_bits;
  int ins_ret = 0;
  if(!kr->restricted || (prefix_i % kr->thread_n) == kr->thread_i){
    ins_ret =  sh_n_add_kmer(kr->sh, source, kmer);
    if(ins_ret > 0) ++(*word_count);
    if(ins_ret == 1) ++(*kmer_count);
  }
  return(ins_ret);
}

void *kmer_reader_read(void *v_args){
  kmer_reader *kr = (kmer_reader*) v_args;
  gzFile gz = gzopen(kr->file_name, "r");
  kseq_t *ks = kseq_init(gz); // remember to gzclose(fp), kseq_destroy(ks)
  size_t read_n = 0;
  int l;
  size_t word_count = 0;
  size_t kmer_count = 0;
  unsigned char *qual;
  unsigned char *seq;
  kmer_iterator km_it;
  kmer_iterator_init(&km_it, kr->k, kr->min_q);
  while( (l = kseq_read(ks)) >= 0 && read_n < kr->max_reads ){
    ++read_n;
    if(l <= kr->k)
      continue;
    qual = (unsigned char*)(ks->qual.l == ks->seq.l ? ks->qual.s : 0);
    seq = (unsigned char*)ks->seq.s;
    uint64_t kmer_f=0;
    uint64_t kmer_r=0;
    // using these iterators makes the code quite a bit slower;
    // I should consider to write a function that procsses one read
    // using the ugly code below... but that's for later.
    if(!kmer_iterator_begin(&km_it, seq, qual, &kmer_f, &kmer_r))
      continue;
    if( kr_add_kmer(kr, kmer_f, kmer_r, &word_count, &kmer_count, kr->source) < 0 )
      continue;
    while(kmer_iterator_next(&km_it, &kmer_f, &kmer_r)){
      if( kr_add_kmer(kr, kmer_f, kmer_r, &word_count, &kmer_count, kr->source) < 0 )
	break;
    }
  }
  gzclose(gz);
  kseq_destroy(ks);
  printf("thread: %d exit.\t%ld\t%ld\n", kr->thread_i, word_count, kmer_count);
  return(v_args);
}

suffix_hash_n* init_kmer_reader_pool(kmer_reader_pool *krp, const char *file_name,
				     int k, uint32_t prefix_bits, size_t max_size, uint32_t thread_n,
				     unsigned char min_q, size_t max_reads,
				     uint32_t source_n, uint32_t source){
  thread_n = (thread_n == 0) ? 1 : thread_n;
  krp->thread_n = thread_n;
  krp->readers = 0;
  uint32_t total_bits = k * 2;
  if(total_bits > SH_KMAX_BITS)
    return(0);
  prefix_bits = (prefix_bits > SH_MAX_PRE_BITS) ? SH_MAX_PRE_BITS : prefix_bits;
  uint32_t suffix_bits = total_bits - prefix_bits;
  // check that suffix bits is not too large; if it is we need to change both suffix and prefix bits
  if(suffix_bits > SH_MAX_SUF_BITS){
    suffix_bits = SH_MAX_SUF_BITS;
    prefix_bits = total_bits - suffix_bits;
  }
  // the suffix hash needs to be kept after all functions have returned
  suffix_hash_n *sh = malloc(sizeof(suffix_hash_n));
  init_suffix_hash_n(sh, source_n, k, prefix_bits, suffix_bits);

  krp->readers = calloc( thread_n, sizeof(kmer_reader) );
  unsigned char restricted = thread_n > 1 ? 1 : 0;
  for(uint32_t i=0; i < thread_n; ++i){
    init_kmer_reader( krp->readers + i, file_name, sh,
		      k, suffix_bits, source, min_q,
		      restricted, thread_n, i, max_reads);
    // start the thread...
    pthread_create( &krp->readers[i].thread, NULL, &kmer_reader_read, &krp->readers[i] );
  }
  return(sh);
}

suffix_hash_n *init_kmer_reader_pool_sh(kmer_reader_pool *krp, const char *file_name,
				      int k, suffix_hash_n *sh, size_t max_size, uint32_t thread_n,
				      unsigned char min_q, size_t max_reads, uint32_t source){
  thread_n = (thread_n == 0) ? 1 : thread_n;
  krp->thread_n = thread_n;
  krp->readers = 0;
  if(k != (sh->suffix_bits + sh->prefix_bits)/2){
    printf("Incompatible arguments: k and total bit numbers do not add up");
    return(sh);
  }
  if(source >= sh->counts_n){
    printf("Value of source is too large\n");
    return(sh);
  }
  uint32_t suffix_bits = sh->suffix_bits;
  krp->readers = calloc( thread_n, sizeof(kmer_reader) );
  unsigned char restricted = thread_n > 1 ? 1 : 0;
  for(uint32_t i=0; i < thread_n; ++i){
    init_kmer_reader( krp->readers + i, file_name, sh, k, suffix_bits, source,
		      min_q, restricted, thread_n, i, max_reads);
    // start the thread...
    pthread_create( &krp->readers[i].thread, NULL, &kmer_reader_read, &krp->readers[i] );
  }
  return(sh);
}

void join_kmer_reader_pool(kmer_reader_pool *krp){
  if(krp->readers == 0)
    return;
  void *ret_value;
  for(uint32_t i=0; i < krp->thread_n; ++i)
    pthread_join(krp->readers[i].thread, &ret_value);
}

// this does not free the suffix_hash; that may be kept
// by an external pointer.
void free_kmer_reader_pool(kmer_reader_pool *krp){
  if(krp->readers == 0)
    return;
  free( krp->readers );
  krp->readers = 0;
}

int seq_kmer_counts(const char* seq, size_t seq_l, int* counts, suffix_hash_n *sh, int k){
  uint64_t off_f = 0;
  uint64_t off_r = 0;
  uint64_t kmer_f = 0;
  uint64_t kmer_r = 0;
  uint64_t kmer = 0;
  uint32_t rc_shift = 64 - (k * 2);
  uint64_t mask = (1ULL << (k * 2)) - 1;
  memset(counts, 0, sizeof(int) * seq_l * sh->counts_n);
  if(k*2 != sh->suffix_bits + sh->prefix_bits)
    return(-1);
  size_t i = 0;
  while(seq[i]){
    if(i == 0 || LC(seq[i]) == 'n'){
      i = init_kmer_qual_2(seq, 0, 0, i, &off_f, &off_r, k);
      kmer_f = off_f & mask;
      kmer_r = off_r >> rc_shift;
      kmer = (kmer_f < kmer_r) ? kmer_f : kmer_r;
      sh_kmer_count_n(sh, kmer, counts + (i-k) * (int)sh->counts_n);
      if(!seq[i])
	break;
    }
    off_f = UPDATE_OFFSET(off_f, seq[i]);
    off_r = UPDATE_OFFSET_RC(off_r, seq[i]);
    kmer_f = off_f & mask;
    kmer_r = off_r >> rc_shift;
    kmer = (kmer_f < kmer_r) ? kmer_f : kmer_r;
    sh_kmer_count_n(sh, kmer, counts + (i-k) * (int)sh->counts_n);
    ++i;
  }
  return(1);
}
  
