#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "kmer_tree.h"


// bits must be below MAX_SUF_BITS
suf_counts *init_suf_counts(uint32_t bits){
  if(bits == 0 || bits > MAX_SUF_BITS)
    return(0);
  suf_counts *sc = calloc(1, sizeof(suf_counts));
  sc->suffix_bits = bits;
  sc->suffix_n = ((size_t)1 << bits);
  sc->suffix_counts = calloc(sc->suffix_n, sizeof(uint32_t));
  return(sc);
}

void free_suf_counts(suf_counts *counts){
  if(!counts)
    return;
  free(counts->suffix_counts);
}

kmer_tree init_kmer_tree(uint32_t prefix_bits, uint32_t suffix_bits, size_t max_size){
  kmer_tree kt;
  kt.suffix_bits = suffix_bits;
  uint32_t total_bits = prefix_bits + suffix_bits;
  kt.kmer_mask = (total_bits >= 64) ? ~(size_t)0 : ((size_t)1 << total_bits) - 1;
  kt.suffix_mask = ((size_t)1 << suffix_bits) - 1;
  kt.prefix_bits = prefix_bits;
  kt.prefixes = 0;
  kt.prefix_n = ((size_t)1 << prefix_bits);
  kt.max_size = max_size;
  kt.allocated = 0;
  kt.max_count = 0;
  kt.max_count_kmer = 0;
  if(sizeof(suf_counts*) * kt.prefix_n <= max_size)
    kt.prefixes = calloc(kt.prefix_n, sizeof(suf_counts*));
  return(kt);
}

void free_kmer_tree(kmer_tree *kt){
  if(!kt->prefixes)
    return;
  for(size_t i=0; i < kt->prefix_n; ++i)
    free_suf_counts(kt->prefixes[i]);
}

// returns negative values for errors
// else, the actual count of the kmer
int add_kmer(kmer_tree *kt, uint64_t kmer){
  kmer &= kt->kmer_mask;
  size_t prefix_i = kmer >> kt->suffix_bits;
  size_t suffix = kmer & kt->suffix_mask;
  if(prefix_i > kt->prefix_n)
    return(-1);
  if(!kt->prefixes[prefix_i]){
    if( (kt->allocated + 1) * (sizeof(suf_counts) + 4 * ((size_t)1 << kt->suffix_bits)) <= kt->max_size )
      kt->prefixes[prefix_i] = init_suf_counts(kt->suffix_bits);
    kt->allocated += kt->prefixes[prefix_i] ? 1 : 0;
  }
  if(!kt->prefixes[prefix_i]){
    printf("allocated prefix: %p : %ld allocated\n", kt->prefixes[prefix_i], kt->allocated);
    printf("prefix / suffix bits %d / %d\n", kt->prefix_bits, kt->suffix_bits);
    printf("estimated memory size: %.2e\n",
	   (double)(kt->allocated + 1) * (sizeof(suf_counts) + 4 * ((size_t)1 << kt->suffix_bits)));
    return(-2);
  }
  uint32_t *count = kt->prefixes[prefix_i]->suffix_counts + suffix;
  ++(*count);
  if(*count > kt->max_count){
    kt->max_count = *count;
    kt->max_count_kmer = kmer;
  }
  return(*count);
}

int kmer_count(kmer_tree *kt, uint64_t kmer){
  kmer &= kt->kmer_mask;
  size_t p_i = kmer >> kt->suffix_bits;
  size_t s_i = kmer & kt->suffix_mask;
  return( kt->prefixes[p_i] == 0 ? 0 : kt->prefixes[p_i]->suffix_counts[s_i] );
}

size_t count_spectrum(kmer_tree *kt, double *counts, uint32_t counts_n){
  uint32_t max_count = counts_n - 1;
  size_t n = 0;
  for(size_t i=0; i < kt->prefix_n; ++i){
    if(!kt->prefixes[i])
      continue;
    for(size_t j=0; j < kt->prefixes[i]->suffix_n; ++j){
      uint32_t count = kt->prefixes[i]->suffix_counts[j];
      counts[ count >= max_count ? max_count : count ]++;
      if(count > 0)
	++n;
    }
  }
  return(n);
}
