#include "kmer_pos.h"
#include "kmer_util.h"
#include "ksort.h"

// For sorting integers.
// #define int_lt(a,b) ((a) < (b))
KSORT_INIT_GENERIC(int)


void clear_kmer_h(khash_t(kmer_h) *hash){
  khiter_t k;
  for(k = kh_begin(hash); k != kh_end(hash); ++k){
    if(kh_exist(hash, k)){
	kv_destroy( kh_val(hash, k).v );
	kh_del( kmer_h, hash, k );
    }
  }
  kh_destroy(kmer_h, hash);
}

void sort_kmer_pos(khash_ptr *hash_ptr){
  khash_t(kmer_h) *hash = hash_ptr->hash;
  khiter_t it;
  for(it = kh_begin(hash); it != kh_end(hash); ++it){
    if(!kh_exist(hash, it))
      continue;
    // consider skipping the exists check;
    // or at least giving a warning message.
    kmer_pos_t kv = kh_val(hash, it);
    // there are several options; but lets try introsort
    ks_introsort(int, kv.v.n, kv.v.a);
  }
}

// Returns 0 if kmer exists in hash; otherwise 1
int kmer_h_insert(uint64_t kmer, int pos, khash_t(kmer_h) *hash){
  int ret =0;
  int new_kmer = 0;
  khiter_t k = kh_get(kmer_h, hash, kmer);
  if(k == kh_end(hash)){
    k = kh_put(kmer_h, hash, kmer, &ret);
    if(k == kh_end(hash)) // error should warn or something
      return(-1);
    kv_init(kh_val(hash, k).v);
    new_kmer = 1;
  }
  kv_push( int, kh_val(hash, k).v, pos );
  kh_val(hash, k).kmer = kmer;
  return(new_kmer);
}

// Returns 0 if kmer does not exit in hash. Otherwise it returns a pointer
// to the kmer_pos_t structure holding the data in the hash. Note that this
// must be used immediately and should not be freed
kmer_pos_t *kmer_pos(khash_t(kmer_h) *hash, uint64_t kmer){
  khiter_t k = kh_get(kmer_h, hash, kmer);
  if(k == kh_end(hash) || !kh_exist(hash, k))
    return(0);
  return( &kh_val(hash, k) );
}

// returns the total number of kmers encountered
// seq: ascii encoded sequence
// k   : the kmer length
// hash: the khash
int seq_to_hash(const char *seq, int k, khash_t(kmer_h) *hash){
  size_t i = 0;
  uint64_t offset = 0;
  int word_count = 0;
  // if k is 32, then (1 << (2*k)) may be undefined behaviour
  // the values that are shifted must be specified as 64 bit;
  // it is possible to do this using 0LL and 1LL, but
  // I'm making it explicit here as I don't think, L, and LL
  // are actually guaranteed to be any particular word length.
  uint64_t one = 1;
  uint64_t zero = 0;
  uint64_t mask = k < 32 ? (one << (2*k)) - 1 : ~zero;
  int ins_ret = 0;
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer(seq, i, &offset, k);
    if(!seq[i])
      break;
    ins_ret = kmer_h_insert( offset & mask, i+1-k, hash );
    if(ins_ret < 0)
      return(ins_ret);
    word_count += ins_ret;
    while(seq[i] && LC(seq[i]) != 'n'){
      offset = UPDATE_OFFSET(offset, seq[i]);
      ++i;
      ins_ret = kmer_h_insert( offset & mask, i+1-k, hash );
      if(ins_ret < 0)
	return(ins_ret);
      word_count += ins_ret;
    }
  }
  return(word_count);
}

// this is a utility function. It simply extends a kvec(int) with alternating values
void pair_positions_push(kmer_pos_t *source, kmer_ppos *dest, int dest_pos){
  if(!source)
    return;
  for(size_t i=0; i < source->v.n; ++i){
    kv_push(int, *dest, dest_pos);
    kv_push(int, *dest, source->v.a[i] );
  }
}

kmer_ppos seq_kmer_positions(khash_t(kmer_h) *hash, const char *seq, int k){
  kmer_ppos pair_positions;
  kv_init(pair_positions);
  // we should consider having an iterator function for going across sequences.
  // but for now I'm just copying the one I have above. I think I might have one
  // somewhere, but, just try things first..
  int i = 0;
  uint64_t offset = 0;
  uint64_t one = 1;
  uint64_t zero = 0;
  uint64_t mask = k < 32 ? (one << (2*k)) - 1 : ~zero;
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer(seq, i, &offset, k);
    if(!seq[i])
      break;
    kmer_pos_t *kpos = kmer_pos(hash, offset & mask);
    pair_positions_push( kpos, &pair_positions, i );
    while(seq[i] && LC(seq[i]) != 'n'){
      offset = UPDATE_OFFSET(offset, seq[i]);
      ++i;
      kmer_pos_t *kpos = kmer_pos(hash, offset & mask);
      pair_positions_push( kpos, &pair_positions, i );
    }
  }
  return(pair_positions);
}
