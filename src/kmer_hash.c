#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <zlib.h>
#include <time.h>
#include "khash.h"
#include "kvec.h"
#include "ksort.h"
#include "kseq.h"
#include "kmer_tree.h"
#include "kmer_reader.h"
#include "suffix_hash.h"
// KSEQ_INIT(gzFile, gzread)

/* // a macro to update an offset as defined here */
/* #define UPDATE_OFFSET(off, c) ( ((off) << 2) | (( (c) >> 1) & 3) ) */
/* #define UPDATE_OFFSET_RC(off, c) ( ((off) >> 2) | ( (( (((uint64_t)(c) >> 1) & 3) + 2) % 4) << 62 )) */
/* #define LC(c) ( (c) | 0x20 ) */
/* #define UC(c) ( (c) & 0xCF ) */
/* #define MAX_K 32 */

#define OPT_KMER 0
#define OPT_POS 1
#define OPT_PAIRS 2
#define OPT_COUNT 3
#define N_OPTS 4
const uint32_t flags[N_OPTS] = {1, 2, 4, 8};
static const char* kmer_pos_fields[N_OPTS] = {"kmer", "pos", "pair.pos", "count"};

// Used to translate offsets to sequences.
const char NUC[4] = {'A', 'C', 'T', 'G'};
const char *kmer_hash_tag = "kmer_hash_250930";
const char *kmer_tree_tag = "kmer_tree_250930";
const char *suffix_hash_tag = "suffix_hash_250930";
const char *suffix_hash_mt_tag = "suffix_hash_mt_250930";
//const char *kmer_count_tag = "kmer_count_250930";

// For sorting integers:
// #define int_lt(a,b) ((a) < (b))
KSORT_INIT_GENERIC(int)


// This will translate into code that defines an
// hash using 64 bit integers as keys
// and the value as a kmer_pos_t
// It is intended for relatively short sequences.
// kvec_t(int)
// the flag value can be set by functions in order to
// mark k-mers that should be treated differently.
// the pos flag may be used to mark positions that
// should be treated differently.
// The pos vector can also be used to hold counts from different
// sources.
typedef struct {
  kvec_t(int) v;
  uint64_t kmer;
  uint64_t kmer_flag;
  //  kvec_t(unsigned char) pos_flag; // usage not yet implemented.
} kmer_pos_t;
KHASH_MAP_INIT_INT64(kmer_h, kmer_pos_t);


void clear_kmer_h(khash_t(kmer_h) *hash){
  khiter_t k;
  for(k = kh_begin(hash); k != kh_end(hash); ++k){
    if(kh_exist(hash, k)){
	kv_destroy( kh_val(hash, k).v );
	//	kv_destroy( kh_val(hash, k).pos_flag );
	kh_del( kmer_h, hash, k );
    }
  } 
}

// The khash will need to be kept as an external pointer
// that points to some form of information:
typedef struct {
  khash_t(kmer_h) *hash;
  int k;
  size_t kmer_count;
  int sorted;
} khash_ptr;

typedef struct {
  kmer_tree tree;
} kmer_tree_ptr;

/// utility function;
SEXP mk_strsxp(const char **words, size_t n){
  SEXP words_r = PROTECT(allocVector(STRSXP, n));
  for(size_t i=0; i < n; ++i)
    SET_STRING_ELT(words_r, i, mkChar(words[i]));
  UNPROTECT(1);
  return(words_r);
}

// a more universal function to extract a pointer
// check that the tag value is correct.
// Returns a null pointer if it encounters an error; does not
// call error as this may be used by the finalise function.
void* extract_ext_ptr(SEXP ptr_r, const char *tag){
  if(TYPEOF(ptr_r) != EXTPTRSXP)
    return(0);
  SEXP tag_r = PROTECT(R_ExternalPtrTag(ptr_r));
  // check the tag;
  if(TYPEOF(tag_r) != STRSXP || length(tag_r) != 1 || strcmp( CHAR(STRING_ELT(tag_r, 0)), tag)){
    UNPROTECT(1);
    return(0);
  }
  UNPROTECT(1);
  return( R_ExternalPtrAddr(ptr_r) );
}

// We need a finalise function to clear resources
// when the external pointer goes out of scope
static void finalise_khash_ptr(SEXP ptr_r){
  void *v_ptr = extract_ext_ptr(ptr_r, kmer_hash_tag);
  //  void *v_ptr = R_ExternalPtrAddr(ptr_r);
  if(!v_ptr)
    return;
  khash_ptr *ptr = (khash_ptr*)v_ptr;
  if(ptr->hash){
    clear_kmer_h(ptr->hash);
    kh_destroy(kmer_h, ptr->hash);
    ptr->hash = 0;
  }
}

static void finalise_kmer_tree_ptr(SEXP ptr_r){
  void *v_ptr = extract_ext_ptr(ptr_r, kmer_tree_tag);
  if(!v_ptr)
    return;
  kmer_tree_ptr *ptr = (kmer_tree_ptr*)v_ptr;
  free_kmer_tree(&ptr->tree);
}

static void finalise_suffix_hash_ptr(SEXP ptr_r){
  void *v_ptr = extract_ext_ptr(ptr_r, suffix_hash_tag);
  if(!v_ptr)
    return;
  suffix_hash *ptr = (suffix_hash*)v_ptr;
  free_suffix_hash(ptr);
  free(ptr);
}

static void finalise_suffix_hash_mt_ptr(SEXP ptr_r){
  void *v_ptr = extract_ext_ptr(ptr_r, suffix_hash_mt_tag);
  if(!v_ptr)
    return;
  suffix_hash_mt *ptr = (suffix_hash_mt*)v_ptr;
  free_suffix_hash_mt(ptr);
  free(ptr);
}


/* size_t skip_n(const char *seq, size_t i){ */
/*   while(seq[i] && LC(seq[i]) == 'n') */
/*     ++i; */
/*   return(i); */
/* } */

/* size_t skip_n_qual(const char *seq, const char *qual, char min_q, size_t i){ */
/*   while(seq[i] && (LC(seq[i]) == 'n' || (qual && qual[i] < min_q))) */
/*     ++i; */
/*   return(i); */
/* } */


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

// sets the seq to the sequence represented by the given
// offset. Assumes that *seq is a char array of size
// k + 1. The function will set seq[k] to 0, though it
// would be more efficient to do so elsewhere.
int kmer_seq(char *seq, unsigned int k, uint64_t offset){
  if(k > MAX_K)
    return(-1);
  size_t mask = 3;
  seq[k] = 0;
  for(ssize_t off=k-1; off >= 0; off--){
    seq[ off ] = NUC[ (offset & mask) ];
    offset >>= 2;
  }
  return(0);
}


// returns the value of j; if it runs out of sequence (i.e. last k-mer)
// or it finds an N. The caller must check the return value.
size_t init_kmer_qual(const char *seq, const char *qual, char min_q, size_t i, unsigned long *offset, int k){
  size_t j = 0;
  while(seq[i]){
    *offset = 0;
    for(j=0; j < k && seq[i+j] && LC(seq[i+j]) != 'n' && (!qual || qual[i+j] >= min_q); ++j){
      *offset = UPDATE_OFFSET(*offset, seq[i+j]);
    }
    if(seq[i+j] == 0 || j == k)
      break;
    // otherwise we hit Ns again;
    i = skip_n_qual(seq, qual, min_q, i + j);
    j = 0;
  }
  return( i + j );
}

// do forward and reverse complement
/* size_t init_kmer_qual_2(const char *seq, const char *qual, char min_q, size_t i, */
/* 			uint64_t *offset, uint64_t *offset_rc, */
/* 			int k){ */
/*   size_t j = 0; */
/*   while(seq[i]){ */
/*     *offset = 0; */
/*     *offset_rc = 0; */
/*     for(j=0; j < k && seq[i+j] && LC(seq[i+j]) != 'n' && (!qual || qual[i+j] >= min_q); ++j){ */
/*       *offset = UPDATE_OFFSET(*offset, seq[i+j]); */
/*       *offset_rc = UPDATE_OFFSET_RC(*offset_rc, seq[i+j]); */
/*     } */
/*     if(seq[i+j] == 0 || j == k) */
/*       break; */
/*     // otherwise we hit Ns again; */
/*     i = skip_n_qual(seq, qual, min_q, i + j); */
/*     j = 0; */
/*   } */
/*   return( i + j ); */
/* } */


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

// The same as kmer_h_insert, but only increments one of several counters
// kmer: a long representing the sequence
// hash: the hash
// source: the source of the kmer (a sequence)
// source_n: the total number of sources
// note that source_n > source must be true
//        
// Returns 1 if kmer is new, otherwise 0
int kmer_count_insert(uint64_t kmer, khash_t(kmer_h) *hash, size_t source, size_t source_n){
  // consider removing the check later if the code works as expected.
  if(source >= source_n){
    warning("source (%ld) equal to or larger than source_n (%ld)", source, source_n);
    return(-1);
  }
  int ret =0;
  int new_kmer = 0;
  khiter_t k = kh_get(kmer_h, hash, kmer);
  if(k == kh_end(hash)){
    k = kh_put(kmer_h, hash, kmer, &ret);
    if(k == kh_end(hash)) // error should warn or something
      return(-1);
    kv_init(kh_val(hash, k).v);
    kh_val(hash, k).kmer = kmer;
    kh_val(hash, k).v.a = malloc(source_n * sizeof(int));
    memset(kh_val(hash, k).v.a, 0, sizeof(int) * source_n);
    kh_val(hash, k).v.m = source_n;
    kh_val(hash, k).v.n = source_n;
    new_kmer = 1;
  }
  kh_val(hash, k).v.a[source]++;
  return(new_kmer);
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

// This is almost identical code to seq_to_hash
// It would be better to refactor, but at the moment I can't
// think of a good function signature that would not be confusing
// seq: ascii encoded sequence
// k   : the kmer length
// hash: the khash
// source: an integer giving the source of the sequence
// source_n: the total number of sequence sources
int seq_to_counts(const char *seq, int k, khash_t(kmer_h) *hash, size_t source, size_t source_n){
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
    ins_ret = kmer_count_insert( offset & mask, hash, source, source_n );
    if(ins_ret < 0)
      return(ins_ret);
    word_count += ins_ret;
    while(seq[i] && LC(seq[i]) != 'n'){
      offset = UPDATE_OFFSET(offset, seq[i]);
      ++i;
      ins_ret = kmer_count_insert( offset & mask, hash, source, source_n );
      if(ins_ret < 0)
	return(ins_ret);
      word_count += ins_ret;
    }
  }
  return(word_count);
}

// like seq_to_counts, but use kmer_tree type of data structure instead
// currently this ignores source_n and source
// as kmer_tree does not support this at the moment
int seq_to_counts_kt(const char *seq, char *qual, char min_q, int k, kmer_tree *kt, size_t source, size_t source_n, size_t *total_counted){
  size_t i = 0;
  uint64_t offset = 0;
  uint64_t offset_rc = 0;
  int word_count = 0;
  // ULL: unsigned long long; necessary to get correct values from shifts
  uint64_t mask = k < 32 ? (1ULL << (2*k)) - 1 : ~0ULL;
  uint32_t rc_shift = 64 - k * 2;
  int ins_ret = 0;
  //  char *kmer_string = calloc(1 + 32, sizeof(char));
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer_qual_2(seq, qual, min_q, i, &offset, &offset_rc, k);
    if(!seq[i])
      break;
    uint64_t off_fwd = offset & mask;
    uint64_t off_rev = (offset_rc >> rc_shift) & mask;
    ins_ret = add_kmer(kt, off_fwd < off_rev ? off_fwd : off_rev);
    if(ins_ret < 0)
      return(ins_ret);
    word_count += (ins_ret == 1 ? 1 : 0);
    while(seq[i] && LC(seq[i]) != 'n' && (!qual || qual[i] > min_q)){
      offset = UPDATE_OFFSET(offset, seq[i]);
      offset_rc = UPDATE_OFFSET_RC(offset_rc, seq[i]);
      off_fwd = mask & offset;
      off_rev = mask & ( offset_rc >> rc_shift );
      ++i;
      ins_ret = add_kmer(kt, off_fwd < off_rev ? off_fwd : off_rev);
      if(ins_ret < 0)
	return(ins_ret);
      ++(*total_counted);
      word_count += (ins_ret == 1 ? 1 : 0);
    }
  }
  //  free(kmer_string);
  return(word_count);
}

// like seq_to_counts_kt, but use suffix_hash type of data structure instead
int seq_to_counts_sh(const char *seq, char *qual, char min_q, int k, suffix_hash *kt, size_t *total_counted){
  size_t i = 0;
  uint64_t offset = 0;
  uint64_t offset_rc = 0;
  int word_count = 0;
  // ULL: unsigned long long; necessary to get correct values from shifts
  uint64_t mask = k < 32 ? (1ULL << (2*k)) - 1 : ~0ULL;
  uint32_t rc_shift = 64 - k * 2;
  int ins_ret = 0;
  //  char *kmer_string = calloc(1 + 32, sizeof(char));
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer_qual_2(seq, qual, min_q, i, &offset, &offset_rc, k);
    if(!seq[i])
      break;
    uint64_t off_fwd = offset & mask;
    uint64_t off_rev = (offset_rc >> rc_shift) & mask;
    ins_ret = sh_add_kmer(kt, off_fwd < off_rev ? off_fwd : off_rev);
    if(ins_ret < 0)
      return(ins_ret);
    word_count += (ins_ret == 1 ? 1 : 0);
    while(seq[i] && LC(seq[i]) != 'n' && (!qual || qual[i] > min_q)){
      offset = UPDATE_OFFSET(offset, seq[i]);
      offset_rc = UPDATE_OFFSET_RC(offset_rc, seq[i]);
      off_fwd = mask & offset;
      off_rev = mask & ( offset_rc >> rc_shift );
      ++i;
      ins_ret = sh_add_kmer(kt, off_fwd < off_rev ? off_fwd : off_rev);
      if(ins_ret < 0)
	return(ins_ret);
      ++(*total_counted);
      word_count += (ins_ret == 1 ? 1 : 0);
    }
  }
  //  free(kmer_string);
  return(word_count);
}

// like seq_to_counts_sh, but use suffix_hash_mt type of data structure instead
int seq_to_counts_sh_mt(const char *seq, char *qual, char min_q, int k, suffix_hash_mt *sh, size_t *total_counted){
  size_t i = 0;
  uint64_t offset = 0;
  uint64_t offset_rc = 0;
  int word_count = 0;
  // ULL: unsigned long long; necessary to get correct values from shifts
  uint64_t mask = k < 32 ? (1ULL << (2*k)) - 1 : ~0ULL;
  uint32_t rc_shift = 64 - k * 2;
  int ins_ret = 0;
  while(seq[i]){
    // pass any potential Ns
    i = init_kmer_qual_2(seq, qual, min_q, i, &offset, &offset_rc, k);
    if(!seq[i])
      break;
    uint64_t off_fwd = offset & mask;
    uint64_t off_rev = (offset_rc >> rc_shift) & mask;
    // this returns 0 on success, 0 otherwise
    ins_ret = suffix_hash_mt_add(sh, off_fwd < off_rev ? off_fwd : off_rev);
    if(ins_ret < 0)
      return(ins_ret);
    //    word_count += (ins_ret == 1 ? 1 : 0);
    while(seq[i] && LC(seq[i]) != 'n' && (!qual || qual[i] > min_q)){
      offset = UPDATE_OFFSET(offset, seq[i]);
      offset_rc = UPDATE_OFFSET_RC(offset_rc, seq[i]);
      off_fwd = mask & offset;
      off_rev = mask & ( offset_rc >> rc_shift );
      ++i;
      ins_ret = suffix_hash_mt_add(sh, off_fwd < off_rev ? off_fwd : off_rev);
      if(ins_ret < 0)
	return(ins_ret);
      ++(*total_counted);
      //      word_count += (ins_ret == 1 ? 1 : 0);
    }
  }
  //  free(kmer_string);
  // this is just 0 if no error.
  return(word_count);
}


// Like seq_to_counts, but minimal quality value within a counted k-mer
// included to reduce number of erroneous calls.
// min_q: taken as is, not converted within this function. 
int seq_to_counts_minq(const char *seq, char *qual, char min_q, int k, khash_t(kmer_h) *hash, size_t source, size_t source_n){
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
    i = init_kmer_qual(seq, qual, min_q, i, &offset, k);
    if(!seq[i])
      break;
    ins_ret = kmer_count_insert( offset & mask, hash, source, source_n );
    if(ins_ret < 0)
      return(ins_ret);
    word_count += ins_ret;
    while(seq[i] && LC(seq[i]) != 'n' && qual[i] >= min_q){
      offset = UPDATE_OFFSET(offset, seq[i]);
      ++i;
      ins_ret = kmer_count_insert( offset & mask, hash, source, source_n );
      if(ins_ret < 0)
	return(ins_ret);
      word_count += ins_ret;
    }
  }
  return(word_count);
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

// COMMENTED OUT: the complexity of task attempted in this
// function requires more concentration than I can give it
// at the moment. It is also likely that it's not really that
// good an idea; I can to some extent just use minimap2 alignments
// or other information.
// But there is something to be said for a non-heuristic function
// for identifying the structure of the underlying sequences.
// construct diagonals from the sequence hash;
// In essence merge diagonals that are consistent;
// extract non central diagonals
/* kvec_t(segment_t) offset_diagonals(const char* seq, khash_ptr *hash_ptr){ */
/*   int k = hash_ptr->k; */
/*   khash_t(kmer_h) *hash = hash_ptr->hash; */
/*   uint64_t offset_i = 0; */
/*   uint64_t offset_j = 0; */
/*   uint64_t mask = k < 32 ? (1 << (2*k)) - 1 : ~0LL; */
/*   kvect_t(segment_t) segments; */
/*   kv_init(segments); */
/*   // then loop through sequence; */
/*   // keep pointers to the previous hash element */
/*   // check if distances are congruent; */
/*   size_t i = 0; */
/*   size_t j = 0; */
/*   //  kvect_t(int) distances; */
/*   //  kvec_init(distances); */
/*   segment_t segment; */
/*   while(seq[i]){ */
/*     i = init_kmer(seq, i, &offset_i, k); */
/*     if(!seq[i]) */
/*       break; */
/*     while(seq[i]){ */
/*       offset_j = offset_i; */
/*       distances.n = 0; */
/*       khiter_t it = kh_get(kmer_h, hash, offset_i); */
/*       if(it != kh_end(hash)){ */
/* 	kvec_t(int) pos = kh_val(hash, it).v; */
/* 	kvec_t(unsigned char) pos_flag = kh_val(hash, it).pos_flag; */
/* 	if(pos_flag.a == 0){ */
/* 	  // we have not come across this kmer before. We will  */
/* 	} */
/* 	segment.x = 1 + i - k; */
/* 	for(size_t m = 1; m < pos.n; ++m){ */
/* 	  // here we should first check whether the position has already */
/* 	  // been incorporated into a diagonal or not. */
/* 	  segment.y = pos.a[m]; */
/* 	  kv_push(segment_t, segments, segment); */
/* 	} */
/*       } */
/*       j = i; */
/*       offset_i = UPDATE_OFFSET(offset_i, seq[i]); */
/*       ++i; */
/*       while(seq[j]){ */
/* 	offset_j = UPDATE_OFFSET(offset_j, seq[j]); */
/* 	khiter_t jt = kh_get(kmer_h, hash, offset); */
/* 	size_t seg_i = 1; */
/* 	int n_merged = 0; */
/* 	if(jt != kh_end(hash)){ */
/* 	  pos = kh_val(hash, jt).v; */
/* 	  for(size_t m = 1; m < pos.n; ++m){ */
/* 	    while(seg_i < segments.n && segments.a[seg_i] < pos.a[m]) */
/* 	      ++seg_i; */
/* 	    if(pos.a[m] - 1 == segments.a[seg_i].y){ */
/* 	      segments.a[seg_i] = pos.a[m]; */
/* 	    } */
/* 	  } */
/* 	} */
/* 	if(n_merged == 0) */
/* 	  break; */
/* 	++j; */
/*       } */
/*     } */
/*   } */
/* } */


khash_ptr* extract_khash_ptr(SEXP ptr_r){
  if(TYPEOF(ptr_r) != EXTPTRSXP)
    error("ptr_r should be an external pointer");
  SEXP tag = PROTECT(R_ExternalPtrTag(ptr_r));
  // check the tag;
  if(TYPEOF(tag) != STRSXP || length(tag) != 1 || strcmp( CHAR(STRING_ELT(tag, 0)), kmer_hash_tag)){
    UNPROTECT(1);
    error("External pointer has incorrect tag");
  }
  khash_ptr *hash_ptr = (khash_ptr*)R_ExternalPtrAddr(ptr_r);
  UNPROTECT(1);
  return(hash_ptr); // which should be checked by the caller
}


SEXP make_kmer_h_index(SEXP seq_r, SEXP k_r, SEXP sort_pos_r){
  if(TYPEOF(seq_r) != STRSXP || length(seq_r) < 1)
    error("seq_r should be a character vector of length at least one");
  if(TYPEOF( k_r ) != INTSXP || length(k_r) < 1 )
    error("k_r must be an integer vector of length at least one");
  if(TYPEOF( sort_pos_r ) != INTSXP || length(k_r) < 1 )
    error("sort_pos_r must be an integer vector of length at least one");
  int k = INTEGER( k_r )[0];
  int sort_pos = asInteger(sort_pos_r);
  if(k < 1 || k > MAX_K)
    error("k must be a positive integer less than 1+MAX_K");
  SEXP seq_r1 = STRING_ELT(seq_r, 0);
  int seq_l = length(seq_r1);
  if(seq_l <= k)
    error("the length of the sequence must be at least k");
  const char *seq = CHAR(seq_r1);
  // Then we need to initialise the khash...
  // calloc initialises memory to 0
  khash_ptr *khash_p = calloc( 1, sizeof(khash_ptr) );
  khash_p->k = k;
  khash_p->hash = kh_init(kmer_h);
  khash_p->kmer_count = seq_to_hash(seq, k, khash_p->hash);
  if(sort_pos)
    sort_kmer_pos(khash_p);
  // and then set up the external pointer:
  SEXP tag = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT( tag, 0, mkChar(kmer_hash_tag));
  // We can consider setting the protect to something else; but
  // I have no use for it presently.
  SEXP prot = R_NilValue;
  SEXP ptr_r = PROTECT(R_MakeExternalPtr(khash_p, tag, prot));
  R_RegisterCFinalizerEx(ptr_r, finalise_khash_ptr, TRUE);
  UNPROTECT(2);
  return(ptr_r);
}

// hash_ptr_r: an external pointer containing a hash_ptr structure
//             if NULL, or not a pointer, then a new hash_ptr will be
//             created.
// params:     k, source, source_n. All integers. 
//             source_n must be positive (>0) and larger than source
// seq_r:      An R STRSXP object. All sequences within it will be added.
SEXP count_kmers(SEXP hash_ptr_r, SEXP params_r, SEXP seq_r){
  if(TYPEOF(seq_r) != STRSXP || length(seq_r) < 1)
    error("seq_r should be a character vector of length at least one");
  if(TYPEOF( params_r ) != INTSXP || length(params_r) != 3 )
    error("k_r must be an integer vector of length 3");
  int *params = INTEGER(params_r);
  int k = params[0];
  int source = params[1];
  int source_n = params[2];
  if(k < 1 || k > MAX_K)
    error("k must be a positive integer less than 1+MAX_K");
  if(source_n < 1 || source >= source_n)
    error("source_n must be larger than 1 and larger than source");
  SEXP ptr_r = hash_ptr_r;
  khash_ptr *khash_p = 0; 
  int protect_n = 0;
  if(ptr_r == R_NilValue){
    khash_p = calloc( 1, sizeof(khash_ptr) );
    khash_p->k = k;
    khash_p->hash = kh_init(kmer_h);
    SEXP tag = PROTECT(mk_strsxp(&kmer_hash_tag, 1));
    ptr_r = PROTECT(R_MakeExternalPtr(khash_p, tag, R_NilValue));
    R_RegisterCFinalizerEx(ptr_r, finalise_khash_ptr, TRUE);
    protect_n = 2;
  }else{
    khash_p = (khash_ptr*)extract_ext_ptr(hash_ptr_r, kmer_hash_tag);
  }
  if(!khash_p)
    error("failed to extract kmer_hash from external pointer");
  if(khash_p->k != k){
    UNPROTECT(protect_n);
    error("mismatch between specified k and that given in the external pointer");
  }
  for(int i=0; i < length(seq_r); ++i){
    const char *seq = CHAR(STRING_ELT(seq_r, i));
    if(length(STRING_ELT(seq_r, i)) <= k)
      continue;
    int kmer_n = seq_to_counts(seq, k, khash_p->hash, (size_t)source, (size_t)source_n);
    if(kmer_n > 0)
      khash_p->kmer_count += kmer_n;
  }
  UNPROTECT(protect_n);
  return(ptr_r);
}


// hash_ptr_r: an external pointer containing a hash_ptr structure
//             if NULL, or not a pointer, then a new hash_ptr will be
//             created.
// params:     k, prefix bits, report_n, max_memory usage, min quality
//             max_read_n
//             All integers; max-memory usage will be considered as
//             gigabytes and multiplied by 2^30
//             min quality will be converted to a phred 33 encoded
//             char; it must not be more than 255-33.
//             max_read_n: the maximal number of reads to process
//             will be coerced to (size_t); -1 is essentially unlimited
//
//             note: dividing up the k-mers into prefix and
//             suffix and using a full table of all suffix
//             counts doesn't actually work as there are too
//             many prefixes. The suffix data structure
//             needs to be sparse as the prefix array of suffixes
//             has not sufficient sparseness to be useful.
//             For now, the best way around this is to not use
//             a very large k, but to 
// fq_file_r:  The name of fastq file; may be compressed.
// Note:   this should work with fasta files as well.
//         KSEQ_INIT(gzFile, gzread)
// at the top of the file: not sure what that does, but I'll check it later.
SEXP count_kmers_fastq(SEXP hash_ptr_r, SEXP params_r, SEXP fq_file_r){
  if(TYPEOF(fq_file_r) != STRSXP || length(fq_file_r) != 1)
    error("fq_file should be a character vector of length at least one");
  if(TYPEOF( params_r ) != INTSXP || length(params_r) != 6 )
    error("k_r must be an integer vector of length 3");
  const char *fq_file = CHAR(STRING_ELT(fq_file_r, 0));
  int *params = INTEGER(params_r);
  int k = params[0];
  size_t report_n = (size_t)params[1];
  if(report_n < 1){
    report_n = 1e6;
    warning("negative or 0 report n value specified: set to 1e6");
  }
  uint32_t prefix_bits = (uint32_t)params[2];
  size_t max_memory = (1ULL << 30) * (size_t)params[3];
  char min_quality = (char)params[4] + '!';
  size_t max_read_n = (size_t)params[5];
  if(k < 1 || k > MAX_K)
    error("k must be a positive integer less than 1+MAX_K");
  uint32_t total_bits = 2 * (uint32_t)k;
  // Note that this is not likely to be optimal, but for testing purposes:
  uint32_t suffix_bits = total_bits - prefix_bits;
  SEXP ptr_r = hash_ptr_r;
  kmer_tree_ptr *ktree_p = 0; 
  int protect_n = 0;
  if(ptr_r == R_NilValue){
    ktree_p = calloc( 1, sizeof(kmer_tree_ptr) );
    ktree_p->tree = init_kmer_tree( prefix_bits, suffix_bits, max_memory );
    SEXP tag = PROTECT(mk_strsxp(&kmer_tree_tag, 1));
    ptr_r = PROTECT(R_MakeExternalPtr(ktree_p, tag, R_NilValue));
    R_RegisterCFinalizerEx(ptr_r, finalise_kmer_tree_ptr, TRUE);
    protect_n = 2;
  }else{
    ktree_p = (kmer_tree_ptr*)extract_ext_ptr(hash_ptr_r, kmer_tree_tag);
  }
  if(!ktree_p)
    error("failed to extract ktree pointer from external pointer");
    
  gzFile fp;
  kseq_t *seq;
  fp = gzopen(fq_file, "r");
  seq = kseq_init(fp);
  int l;
  size_t n_reads = 0;
  // print timing information every report_n reads
  // include the most frequent kmer
  char *kmer_string = calloc(k+1, sizeof(char));
  size_t last_kmer_count = 0;
  size_t total_kmer_count = 0;
  size_t last_allocated_n = 0;
  clock_t read_beg, read_ticks=0;
  //  read_ticks = 0;
  clock_t time_end = 0;
  clock_t time_beg = read_beg = clock();
  size_t total_counted = 0;
  while( (l = kseq_read(seq)) >= 0 && n_reads < max_read_n){
    read_ticks += (clock() - read_beg);
    ++n_reads;
    if(seq->seq.l <= k)
      continue;
    int kmer_n = seq_to_counts_kt(seq->seq.s, seq->qual.s, min_quality, k, &ktree_p->tree, 0, 1, &total_counted);
    if(kmer_n < 0){
      Rprintf("received error code: %d at read: %ld\n", kmer_n, n_reads);
      break;
    }
    total_kmer_count += kmer_n;
    if(n_reads % report_n == 0){
      time_end = clock();
      // Note: CLOCKS_PER_SEC is not the same as the CPU frequency
      // instead it is set to 1e6
      clock_t clicks = time_end - time_beg;
      Rprintf("%ld reads -> %.2e kmers (%.2e new): %.3e clicks %ld sec / %ld reads:  %.3e clicks / read (IO: %.3e)\ntotal allocated: %ld (%.2e)  new: %ld\n",
	      n_reads, (double)total_kmer_count, (double)total_kmer_count - (double)last_kmer_count,
	      (double)clicks, clicks / CLOCKS_PER_SEC, report_n, (double)clicks / (double)report_n, (double)read_ticks,
	      ktree_p->tree.allocated, (double)total_kmer_count / (double)ktree_p->tree.allocated,
	      ktree_p->tree.allocated - last_allocated_n);
      last_allocated_n = ktree_p->tree.allocated;
      Rprintf("prefix / suffix bits: %d / %d\n", ktree_p->tree.prefix_bits, ktree_p->tree.suffix_bits);
      Rprintf("Estimated memory use: %.2e\n", (double)( ktree_p->tree.allocated * 4 * ((size_t)1 << ktree_p->tree.suffix_bits)));
      Rprintf("Total number of words counted %ld max count: %d\n", total_counted, ktree_p->tree.max_count);
      kmer_seq( kmer_string, k, ktree_p->tree.max_count_kmer );
      Rprintf("Most common k-mer: %s\n\n", kmer_string);
      read_ticks = 0;
      last_kmer_count = total_kmer_count;
      time_beg = clock();
    }
    read_beg = clock();
  }
  free(kmer_string);
  kseq_destroy(seq);
  gzclose(fp);
  UNPROTECT(protect_n);
  return(ptr_r);
}

// like count_kmers_fastq
// but uses an array of hashes to hold the k-mers...
SEXP count_kmers_fastq_sh(SEXP hash_ptr_r, SEXP params_r, SEXP fq_file_r){
  if(TYPEOF(fq_file_r) != STRSXP || length(fq_file_r) != 1)
    error("fq_file should be a character vector of length at least one");
  if(TYPEOF( params_r ) != INTSXP || length(params_r) != 6 )
    error("k_r must be an integer vector of length 3");
  const char *fq_file = CHAR(STRING_ELT(fq_file_r, 0));
  int *params = INTEGER(params_r);
  int k = params[0];
  size_t report_n = (size_t)params[1];
  if(report_n < 1){
    report_n = 1e6;
    warning("negative or 0 report n value specified: set to 1e6");
  }
  uint32_t prefix_bits = (uint32_t)params[2];
  size_t max_memory = (1ULL << 30) * (size_t)params[3];
  char min_quality = (char)params[4] + '!';
  size_t max_read_n = (size_t)params[5];
  if(k < 1 || k > MAX_K)
    error("k must be a positive integer less than 1+MAX_K");
  uint32_t total_bits = 2 * (uint32_t)k;
  uint32_t suffix_bits = total_bits - prefix_bits;
  SEXP ptr_r = hash_ptr_r;
  int protect_n = 0;
  suffix_hash *suf_hash = 0; 
  if(ptr_r == R_NilValue){
    suf_hash = calloc(1, sizeof(suffix_hash));
    *suf_hash = init_suffix_hash(prefix_bits, suffix_bits, max_memory);
    SEXP tag = PROTECT(mk_strsxp(&suffix_hash_tag, 1));
    ptr_r = PROTECT(R_MakeExternalPtr(suf_hash, tag, R_NilValue));
    R_RegisterCFinalizerEx(ptr_r, finalise_suffix_hash_ptr, TRUE);
    protect_n = 2;
  }else{
    suf_hash = (suffix_hash*)extract_ext_ptr(ptr_r, suffix_hash_tag);
  }
  if(!suf_hash)
    error("Unable to extract suffix_hash from external pointer");
  gzFile fp;
  kseq_t *seq;
  fp = gzopen(fq_file, "r");
  seq = kseq_init(fp);
  int l;
  size_t n_reads = 0;
  // print timing information every report_n reads
  // include the most frequent kmer
  char *kmer_string = calloc(k+1, sizeof(char));
  size_t last_kmer_count = 0;
  size_t total_kmer_count = 0;
  size_t last_allocated_n = 0;
  clock_t read_beg, read_ticks=0;
  //  read_ticks = 0;
  clock_t time_end = 0;
  clock_t time_beg = read_beg = clock();
  size_t total_counted = 0;
  while( (l = kseq_read(seq)) >= 0 && n_reads < max_read_n){
    read_ticks += (clock() - read_beg);
    ++n_reads;
    if(seq->seq.l <= k)
      continue;
    int kmer_n = seq_to_counts_sh(seq->seq.s, seq->qual.s, min_quality, k, suf_hash, &total_counted);
    if(kmer_n < 0){
      Rprintf("received error code: %d at read: %ld\n", kmer_n, n_reads);
      break;
    }
    total_kmer_count += kmer_n;
    if(n_reads % report_n == 0){
      time_end = clock();
      // Note: CLOCKS_PER_SEC is not the same as the CPU frequency
      // instead it is set to 1e6
      clock_t clicks = time_end - time_beg;
      Rprintf("%ld reads -> %.2e kmers (%.2e new): %.3e clicks %ld sec / %ld reads:  %.3e clicks / read (IO: %.3e)\ntotal allocated: %ld (%.2e)  new: %ld\n",
	      n_reads, (double)total_kmer_count, (double)total_kmer_count - (double)last_kmer_count,
	      (double)clicks, clicks / CLOCKS_PER_SEC, report_n, (double)clicks / (double)report_n, (double)read_ticks,
	      suf_hash->allocated, (double)total_kmer_count / (double)suf_hash->allocated,
	      suf_hash->allocated - last_allocated_n);
      last_allocated_n = suf_hash->allocated;
      Rprintf("prefix / suffix bits: %d / %d\n", suf_hash->prefix_bits, suf_hash->suffix_bits);
      Rprintf("Estimated memory use: %.2e\n", (double)( suf_hash->allocated * 4 * ((size_t)1 << suf_hash->suffix_bits)));
      Rprintf("Total number of words counted %ld max count: %d\n", total_counted, suf_hash->max_count);
      kmer_seq( kmer_string, k, suf_hash->max_count_kmer );
      Rprintf("Most common k-mer: %s\n\n", kmer_string);
      read_ticks = 0;
      last_kmer_count = total_kmer_count;
      time_beg = clock();
    }
    read_beg = clock();
  }
  free(kmer_string);
  kseq_destroy(seq);
  gzclose(fp);
  UNPROTECT(protect_n);
  return(ptr_r);
}

// use a reader pool
// params should be k, prefix_bits, min_qual, thread_n, max_read_n, max_memory
// Currently we do not use the hash_ptr_r; that should be done later
SEXP count_kmers_fastq_sh_rp(SEXP hash_ptr_r, SEXP params_r, SEXP fq_file_r){
  if(TYPEOF(fq_file_r) != STRSXP || length(fq_file_r) != 1)
    error("fq_file should be a character vector of length at least one");
  if(TYPEOF( params_r ) != INTSXP || length(params_r) != 6 )
    error("k_r must be an integer vector of length 6 (k, prefix_bits, min_q, thread_n, max_reads, max_mem");
  const char *fq_file = CHAR(STRING_ELT(fq_file_r, 0));
  int *params = INTEGER(params_r);
  int k = params[0];
  uint32_t prefix_bits = (uint32_t)params[1];
  unsigned char min_q = '!' + (unsigned char)params[2];
  uint32_t thread_n = (uint32_t)params[3];
  size_t max_read_no = (size_t)params[4];
  size_t max_memory = (1ULL << 30) * (size_t)params[5];
  if(k < 1 || k > MAX_K)
    error("k must be a positive integer less than 1+MAX_K");
  // this should run everything in one go
  kmer_reader_pool *krp = malloc(sizeof(kmer_reader_pool));
  // check if we have a valid pointer or not:
  suffix_hash *sh = extract_ext_ptr(hash_ptr_r, suffix_hash_tag);
  int sh_is_new = (sh == 0) ? 1 : 0;
  if(!sh)
    sh = init_kmer_reader_pool(krp, fq_file, k, prefix_bits, max_memory, thread_n, min_q, max_read_no);
  else
    sh = init_kmer_reader_pool_sh(krp, fq_file, k, sh, max_memory, thread_n, min_q, max_read_no);
  Rprintf("read threads started: returned suffix_hash: %p\n", sh);
  if(!sh)
    error("Did not obtain a suffix_hash");
  join_kmer_reader_pool(krp);
  // we could obtain some information from the read pool at this point. But lets not bother.
  free_kmer_reader_pool(krp);
  free(krp);
  if(!sh_is_new)
    return(hash_ptr_r);
  // otherwise make a new external pointer.
  SEXP tag = PROTECT(mk_strsxp(&suffix_hash_tag, 1));
  SEXP ptr_r = PROTECT(R_MakeExternalPtr(sh, tag, R_NilValue));
  R_RegisterCFinalizerEx(ptr_r, finalise_suffix_hash_ptr, TRUE);
  UNPROTECT(2);
  return(ptr_r);
}

SEXP seq_kmer_depth_sh(SEXP hash_ptr_r, SEXP seq_r, SEXP k_r){
  suffix_hash *sh = extract_ext_ptr(hash_ptr_r, suffix_hash_tag);
  if(!sh)
    error("unable to obtain suffix_hash from external pointer");
  if(TYPEOF(k_r) != INTSXP || length(k_r) != 1)
    error("k_r should be a single integer");
  int k = asInteger(k_r);
  if(TYPEOF(seq_r) != STRSXP || length(seq_r) != 1)
    error("seq_r should be a character vector of length 1");
  SEXP seq_rr = STRING_ELT(seq_r, 0);
  size_t seq_l = (size_t)length(seq_rr);
  const char *seq = CHAR(seq_rr);
  SEXP counts_r = PROTECT(allocVector(INTSXP, seq_l));
  int *counts = INTEGER(counts_r);
  if( seq_kmer_counts(seq, seq_l, counts, sh, k) != 1 ){
    UNPROTECT(1);
    error("Receieved error from seq_kmer_counts");
  }
  UNPROTECT(1);
  return(counts_r);
}
  


// A copy of count_kmers_fastq_sh
// These functions should be restructured to avoid the current mess
// params: k, report_n, prefix_bits, thread_n, min_quality, max_read_n
//         queue_buffer_size
SEXP count_kmers_fastq_sh_mt(SEXP hash_ptr_r, SEXP params_r, SEXP fq_file_r){
  if(TYPEOF(fq_file_r) != STRSXP || length(fq_file_r) != 1)
    error("fq_file should be a character vector of length at least one");
  if(TYPEOF( params_r ) != INTSXP || length(params_r) != 7 )
    error("k_r must be an integer vector of length 7");
  const char *fq_file = CHAR(STRING_ELT(fq_file_r, 0));
  int *params = INTEGER(params_r);
  int k = params[0];
  size_t report_n = (size_t)params[1];
  if(report_n < 1){
    report_n = 1e6;
    warning("negative or 0 report n value specified: set to 1e6");
  }
  uint32_t prefix_bits = (uint32_t)params[2];
  //  size_t max_memory = (1ULL << 30) * (size_t)params[3];
  size_t thread_n = params[3];
  char min_quality = (char)params[4] + '!';
  size_t max_read_n = (size_t)params[5];
  if(k < 1 || k > MAX_K)
    error("k must be a positive integer less than 1+MAX_K");
  // the limit of the queue_buffer size should depend on the numbe of threads;
  // if a small number of threads then a large buffer size is OK. But,
  // to some extent the buffer size doesn't really matter as we will almost certainly
  // be filling them up, and once it's full it will just keep filling up.
  size_t q_buffer_size = (size_t)params[6];
  if(q_buffer_size < 1 || q_buffer_size > 1e6) 
    error("q_buffer size too small or too large");
  uint32_t total_bits = 2 * (uint32_t)k;
  uint32_t suffix_bits = total_bits - prefix_bits;
  // These should be checked here..
  if(total_bits > SH_KMAX_BITS || prefix_bits > SH_MAX_PRE_BITS || suffix_bits > SH_MAX_SUF_BITS)
    error("too many bits somewhere");
  SEXP ptr_r = hash_ptr_r;
  int protect_n = 0;
  suffix_hash_mt *suf_hash = 0; 
  if(ptr_r == R_NilValue){
    suf_hash = calloc(1, sizeof(suffix_hash_mt));
    init_suffix_hash_mt(suf_hash, prefix_bits, suffix_bits, thread_n, q_buffer_size);
    SEXP tag = PROTECT(mk_strsxp(&suffix_hash_mt_tag, 1));
    ptr_r = PROTECT(R_MakeExternalPtr(suf_hash, tag, R_NilValue));
    R_RegisterCFinalizerEx(ptr_r, finalise_suffix_hash_mt_ptr, TRUE);
    protect_n = 2;
  }else{
    suf_hash = (suffix_hash_mt*)extract_ext_ptr(ptr_r, suffix_hash_mt_tag);
  }
  if(!suf_hash)
    error("Unable to extract suffix_hash from external pointer");
  gzFile fp;
  kseq_t *seq;
  fp = gzopen(fq_file, "r");
  seq = kseq_init(fp);
  int l;
  size_t n_reads = 0;
  // print timing information every report_n reads
  // include the most frequent kmer
  //  char *kmer_string = calloc(k+1, sizeof(char));
  clock_t time_end = 0;
  clock_t time_beg = clock();
  size_t total_counted = 0;
  while( (l = kseq_read(seq)) >= 0 && n_reads < max_read_n){
    ++n_reads;
    if(seq->seq.l <= k)
      continue;
    int ret = seq_to_counts_sh_mt(seq->seq.s, seq->qual.s, min_quality, k, suf_hash, &total_counted);
    if(ret < 0){
      Rprintf("received error code: %d at read: %ld\n", ret, n_reads);
      break;
    }
    if(n_reads % report_n == 0){
      time_end = clock();
      // Note: CLOCKS_PER_SEC is not the same as the CPU frequency
      // instead it is set to 1e6
      clock_t clicks = time_end - time_beg;
      Rprintf("%ld reads: %ld reads in %.3e clicks (%f seconds): %.3e reads / second\n",
	      n_reads, report_n, (double)clicks, (double)clicks / (double)CLOCKS_PER_SEC,
	      (double)report_n / ((double)clicks / (double)CLOCKS_PER_SEC));
      Rprintf("Total of %ld words counted\n", total_counted);
      time_beg = clock();
    }
  }
  suffix_hash_mt_join(suf_hash);
  kseq_destroy(seq);
  gzclose(fp);
  UNPROTECT(protect_n);
  return(ptr_r);
}


SEXP kmer_spectrum_ktree(SEXP ext_ptr, SEXP max_count_r){
  kmer_tree_ptr *ktree_p = (kmer_tree_ptr*)extract_ext_ptr(ext_ptr, kmer_tree_tag);
  if(ktree_p == 0)
    error("unable to obtain a kmer_tree_ptr from external pointer");
  kmer_tree kt = ktree_p->tree;
  if(TYPEOF(max_count_r) != INTSXP || length(max_count_r) != 1)
    error("max_count_r should be a single integer");
  int max_count = asInteger(max_count_r);
  if(max_count < 1 || max_count > (1 << 30))
    error("Unsuitable value of max_count");
  SEXP counts_r = PROTECT(allocVector(REALSXP, max_count + 1));
  double *counts = REAL(counts_r);
  memset(counts, 0, sizeof(double) * (max_count + 1));
  count_spectrum(&kt, counts, max_count + 1);
  UNPROTECT(1);
  return(counts_r);
}

SEXP kmer_spectrum_suffix_hash(SEXP ext_ptr, SEXP max_count_r){
  suffix_hash *suf_hash = (suffix_hash*)extract_ext_ptr(ext_ptr, suffix_hash_tag);
  if(suf_hash == 0)
    error("unable to obtain a suffix_hash from external pointer");
  if(TYPEOF(max_count_r) != INTSXP || length(max_count_r) != 1)
    error("max_count_r should be a single integer");
  int max_count = asInteger(max_count_r);
  if(max_count < 1 || max_count > (1 << 30))
    error("Unsuitable value of max_count");
  SEXP counts_r = PROTECT(allocVector(REALSXP, max_count + 1));
  double *counts = REAL(counts_r);
  memset(counts, 0, sizeof(double) * (max_count + 1));
  sh_count_spectrum(suf_hash, counts, max_count + 1);
  UNPROTECT(1);
  return(counts_r);
}


// To count kmers from fastq files I can use kseq.h
// this allows me to read one sequence at a time from compressed fastq / fasta files
// in theory, this means that I could use the quality information present in the
// read data.
// See example usage of kseq.h at
// https://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser


// Return list containing optionally:
// character vector: giving the kmers
// a matrix giving the index of the kmer followed
// by the position; kmer_i, pos, kmer_i pos
// This is not the most efficient of data structures, but it will
// suffice for now. We can simply use a kvec_h(int)
SEXP kmer_positions(SEXP ptr_r, SEXP opt_flag_r){
  khash_ptr *h_ptr = extract_khash_ptr(ptr_r);
  // we have to trust this pointer; We could consider making it's
  // first element a magic number, but not much point as that
  // is what we try to do in the tag.
  if(TYPEOF(opt_flag_r) != INTSXP || length(opt_flag_r) != 1)
    error("opt_flag_r should be an integer vector of length 1");
  uint32_t opt_flag = asInteger(opt_flag_r);
  // returning the kmer sequences is expensive;
  SEXP ret_data = PROTECT(allocVector(VECSXP, N_OPTS));
  setAttrib( ret_data, R_NamesSymbol, mk_strsxp(kmer_pos_fields, N_OPTS) );

  // we define both kmer_pos and kmer_pair_pos
  // even if we are not going to fill them. 
  // to reduce the number of conditionals

  //kmer_pos: a matrix with two rows: i, pos
  kvec_t(int) kmer_pos;
  kv_init(kmer_pos);

  //kmer_pair_pos: a matrix with three rows: i, pos_1, pos_2
  // only defined where pos_1 and pos_2 are different
  kvec_t(int) kmer_pair_pos;
  kv_init(kmer_pair_pos);

  // the number of occurences of each kmer; to obtain these in `R` would
  // require a tapply across the positions table; which can be slow.
  kvec_t(int) kmer_counts;
  kv_init(kmer_counts);
  
  khiter_t k_it;
  khash_t(kmer_h) *hash = h_ptr->hash;
  int k=h_ptr->k;
  char *buffer = malloc(k + 1);
  buffer[k] = 0;
  size_t hash_n = kh_size(hash);
  SEXP kmers_r = R_NilValue;
  if(opt_flag & flags[OPT_KMER]){
    SET_VECTOR_ELT( ret_data, 0, allocVector(STRSXP, hash_n));
    kmers_r = VECTOR_ELT(ret_data, 0);
  }
  int i=0;
  for(k_it = kh_begin(hash); k_it != kh_end(hash); ++k_it){
    if(kh_exist(hash, k_it)){
      //      uint64_t key = kh_key(hash, k_it);
      kmer_pos_t kv = kh_val(hash, k_it);
      //      kmer_seq(buffer, k, key);
      if(opt_flag & flags[OPT_KMER]){
	kmer_seq(buffer, k, kv.kmer);
	SET_STRING_ELT(kmers_r, i, mkChar(buffer));
      }
      if(opt_flag & flags[OPT_COUNT])
	kv_push(int, kmer_counts, kv.v.n);
      ++i;
      if((opt_flag & (flags[OPT_POS] | flags[OPT_PAIRS])) == 0)
	continue;
      for(int j=0; j < kv.v.n; ++j){
	if(opt_flag & flags[OPT_POS]){
	  kv_push(int, kmer_pos, i);
	  kv_push(int, kmer_pos, kv.v.a[j]); // kv_A(kv.v, j));
	}
	if(opt_flag & flags[OPT_PAIRS]){
	  for(int k=j+1; k < kv.v.n; ++k){
	    // It would be better to do this with a single push of a struct
	    // I should change the data structures to do this.
	    kv_push(int, kmer_pair_pos, i);
	    kv_push(int, kmer_pair_pos, kv.v.a[j]);
	    kv_push(int, kmer_pair_pos, kv.v.a[k]);
	  }
	}
      }
    }
  }
  int *pos = 0;
  int *pair_pos = 0;
  if(opt_flag & flags[OPT_POS]){
    SET_VECTOR_ELT(ret_data, OPT_POS, allocMatrix(INTSXP, 2, kmer_pos.n / 2));
    pos = INTEGER(VECTOR_ELT(ret_data, OPT_POS));
    memcpy( pos, kmer_pos.a, kmer_pos.n * sizeof(int));
  }
  if(opt_flag & flags[OPT_PAIRS]){
    SET_VECTOR_ELT(ret_data, OPT_PAIRS, allocMatrix(INTSXP, 3, kmer_pair_pos.n /3));
    pair_pos = INTEGER(VECTOR_ELT(ret_data, OPT_PAIRS));
    memcpy(pair_pos, kmer_pair_pos.a, kmer_pair_pos.n * sizeof(int));
  }
  if(opt_flag & flags[OPT_COUNT]){
    SET_VECTOR_ELT(ret_data, OPT_COUNT, allocVector(INTSXP, kmer_counts.n));
    memcpy(INTEGER(VECTOR_ELT(ret_data, OPT_COUNT)), kmer_counts.a, kmer_counts.n * sizeof(int));
  }
  free(buffer);
  kv_destroy(kmer_pos);
  kv_destroy(kmer_pair_pos);
  kv_destroy(kmer_counts);
  UNPROTECT(1);
  return(ret_data);
}

SEXP kmer_pair_pos(SEXP ptr_a, SEXP ptr_b){
  khash_ptr *a = extract_khash_ptr(ptr_a);
  khash_ptr *b = extract_khash_ptr(ptr_b);
  // collect points in kvec_int
  kvec_t(int) pairs;
  kv_init(pairs);
  khiter_t it_a;
  khiter_t it_b;
  for(it_a = kh_begin(a->hash); it_a != kh_end(a->hash); ++it_a){
    Rprintf( "%ld\n", kh_key(a->hash, it_a));
    it_b = kh_get(kmer_h, b->hash, kh_key(a->hash, it_a));
    if(kh_exist(b->hash, it_b)){
      kmer_pos_t av = kh_val(a->hash, it_a);
      kmer_pos_t bv = kh_val(b->hash, it_b);
      Rprintf("a: %ld  b: %ld\n", av.v.n, bv.v.n);
      for(size_t i=0; i < av.v.n; ++i){
	for(size_t j=0; j < bv.v.n; ++j){
	  kv_push(int, pairs, av.v.a[i]);
	  kv_push(int, pairs, bv.v.a[j]);
	}
      }
    }
  }
  SEXP ret_ptr = PROTECT(allocMatrix(INTSXP, 2, pairs.n/2));
  int *ret_data = INTEGER(ret_ptr);
  memcpy(ret_data, pairs.a, pairs.n * sizeof(int));
  kv_destroy(pairs);
  UNPROTECT(1);
  return(ret_ptr);
}

static const R_CallMethodDef callMethods[] = {
	      {"make_kmer_h_index", (DL_FUNC)&make_kmer_h_index, 3},
	      {"count_kmers", (DL_FUNC)&count_kmers, 3},
	      {"count_kmers_fastq", (DL_FUNC)&count_kmers_fastq, 3},
	      {"count_kmers_fastq_sh", (DL_FUNC)&count_kmers_fastq_sh, 3},
	      {"count_kmers_fastq_sh_mt", (DL_FUNC)&count_kmers_fastq_sh_mt, 3},
	      {"count_kmers_fastq_sh_rp", (DL_FUNC)&count_kmers_fastq_sh_rp, 3},
	      {"seq_kmer_depth_sh", (DL_FUNC)&seq_kmer_depth_sh, 3},
	      {"kmer_spectrum_ktree", (DL_FUNC)&kmer_spectrum_ktree, 2},
	      {"kmer_spectrum_suffix_hash", (DL_FUNC)&kmer_spectrum_suffix_hash, 2},
	      {"kmer_positions", (DL_FUNC)&kmer_positions, 2},
	      {"kmer_pair_pos", (DL_FUNC)&kmer_pair_pos, 2},
	      {NULL, NULL, 0}
};

void R_init_kmer_hash(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
