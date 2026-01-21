#include <stdio.h>
#include <zlib.h>
#include "kmer_reader.h"

// KSEQ_INIT(gzFile, gzread)

size_t skip_n_qual(const char *seq, const char *qual, char min_q, size_t i){
  while(seq[i] && (LC(seq[i]) == 'n' || (qual && qual[i] < min_q)))
    ++i;
  return(i);
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
    if(seq[i+j] == 0 || j == k)
      break;
    // otherwise we hit Ns again;
    i = skip_n_qual(seq, qual, min_q, i + j);
    j = 0;
  }
  return( i + j );
}


void init_kmer_reader(kmer_reader *kr, const char *file, suffix_hash *sh, int k, uint32_t suffix_bits,
		      unsigned char min_q, unsigned char restricted, uint32_t thread_n,
		      uint32_t thread_i, size_t max_reads){
  kr->file_name = file;
  kr->sh = sh;
  kr->k = k;
  kr->suffix_bits = suffix_bits;
  kr->min_q = min_q;
  kr->restricted = restricted;
  kr->thread_n = thread_n;
  kr->thread_i = thread_i;
  kr->max_reads = max_reads;
  // the pthread_t thread is initialised when pthread_create is called
}

void *kmer_reader_read(void *v_args){
  kmer_reader *kr = (kmer_reader*) v_args;
  gzFile gz = gzopen(kr->file_name, "r");
  kseq_t *ks = kseq_init(gz); // remember to gzclose(fp), kseq_destroy(ks)
  size_t read_n = 0;
  int l;
  uint64_t off_fwd = 0;
  uint64_t off_rev = 0;
  uint64_t offset = 0;
  uint64_t offset_rc = 0;
  uint64_t kmer = 0;
  uint64_t mask = (1ULL << (kr->k*2)) - 1;
  size_t word_count = 0;
  size_t kmer_count = 0;
  char *qual;
  char *seq;
  uint32_t rc_shift = 64 - kr->k * 2;
  int ins_ret = 0;
  uint64_t prefix_i = 0;
  while( (l = kseq_read(ks)) >= 0 && read_n < kr->max_reads ){
    ++read_n;
    if(l <= kr->k)
      continue;
    // The following would be better off as a separate function:
    // seq_to_kmer( struct seq_k )
    // where seq_k, contains the kseq pointer as well as the state (i.e. position
    // and if initialised)
    size_t i = 0;
    qual = ks->qual.l == ks->seq.l ? ks->qual.s : 0;
    seq = ks->seq.s;
    while(seq[i]){
    // pass any potential Ns and low quality regions
      i = init_kmer_qual_2(seq, qual, kr->min_q, i, &offset, &offset_rc, kr->k);
      if(!seq[i])
	break;
      off_fwd = offset & mask;
      off_rev = (offset_rc >> rc_shift) & mask;
      kmer = off_fwd < off_rev ? off_fwd : off_rev;
      prefix_i = kmer >> kr->suffix_bits;
      if(!kr->restricted || (prefix_i % kr->thread_n) == kr->thread_i){
	ins_ret = sh_add_kmer(kr->sh, kmer);
	if(ins_ret < 0)
	  break;
	++word_count;
	if(ins_ret == 1)
	  ++kmer_count;
      }
      while(seq[i] && LC(seq[i]) != 'n' && (!qual || qual[i] > kr->min_q)){
	offset = UPDATE_OFFSET(offset, seq[i]);
	offset_rc = UPDATE_OFFSET_RC(offset_rc, seq[i]);
	++i;
	off_fwd = mask & offset;
	off_rev = mask & ( offset_rc >> rc_shift );
	kmer = off_fwd < off_rev ? off_fwd : off_rev;
	prefix_i = kmer >> kr->suffix_bits;
	if(!kr->restricted || (prefix_i % kr->thread_n) == kr->thread_i){
	  ins_ret = sh_add_kmer(kr->sh, kmer);
	  if(ins_ret < 0)
	    break; // this will only break out of the current region; not the whole sequence.
	  ++word_count;
	  if(ins_ret == 1)
	    ++kmer_count;
	}
      }
    }
  }
  gzclose(gz);
  kseq_destroy(ks);
  printf("thread: %d exit.\t%ld\t%ld\n", kr->thread_i, word_count, kmer_count);
  return(v_args);
}

suffix_hash *init_kmer_reader_pool(kmer_reader_pool *krp, const char *file_name,
				   int k, uint32_t prefix_bits, size_t max_size, uint32_t thread_n,
				   unsigned char min_q, size_t max_reads){
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
  suffix_hash *sh = malloc(sizeof(suffix_hash));
  init_suffix_hash_p(sh, prefix_bits, suffix_bits, max_size);
  krp->readers = calloc( thread_n, sizeof(kmer_reader) );
  unsigned char restricted = thread_n > 1 ? 1 : 0;
  for(uint32_t i=0; i < thread_n; ++i){
    init_kmer_reader( krp->readers + i, file_name, sh, k, suffix_bits, min_q, restricted, thread_n, i, max_reads);
    // start the thread...
    pthread_create( &krp->readers[i].thread, NULL, &kmer_reader_read, &krp->readers[i] );
  }
  return(sh);
}

suffix_hash *init_kmer_reader_pool_sh(kmer_reader_pool *krp, const char *file_name,
				      int k, suffix_hash *sh, size_t max_size, uint32_t thread_n,
				      unsigned char min_q, size_t max_reads){
  thread_n = (thread_n == 0) ? 1 : thread_n;
  krp->thread_n = thread_n;
  krp->readers = 0;
  if(k != (sh->suffix_bits + sh->prefix_bits)/2){
    printf("Incompatible arguments: k and total bit numbers do not add up");
    return(sh);
  }
  uint32_t prefix_bits = sh->prefix_bits;
  uint32_t suffix_bits = sh->suffix_bits;
  krp->readers = calloc( thread_n, sizeof(kmer_reader) );
  unsigned char restricted = thread_n > 1 ? 1 : 0;
  for(uint32_t i=0; i < thread_n; ++i){
    init_kmer_reader( krp->readers + i, file_name, sh, k, suffix_bits, min_q, restricted, thread_n, i, max_reads);
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
}

  
