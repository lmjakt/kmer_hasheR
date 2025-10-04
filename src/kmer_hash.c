#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include "khash.h"
#include "kvec.h"
#include "ksort.h"

// a macro to update an offset as defined here
#define UPDATE_OFFSET(off, c) ( ((off) << 2) | (( (c) >> 1) & 3) )
#define LC(c) ( (c) | 0x20 )
#define UC(c) ( (c) & 0xCF )
#define MAX_K 32

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

// For sorting integers:
// #define int_lt(a,b) ((a) < (b))
KSORT_INIT_GENERIC(int)


// This will translate into code that defines an
// hash using 64 bit integers as keys and
// kvec_t(int) as the values...
// the flag value can be set by functions in order to
// mark k-mers that should be treated differently.
// the pos flag may be used to mark positions that
// should be treated differently.
typedef struct {
  kvec_t(int) v;
  uint64_t kmer;
  uint64_t kmer_flag;
  //  kvec_t(unsigned char) pos_flag; // usage not yet implemented.
} kmer_pos_t;
KHASH_MAP_INIT_INT64(kmer_h, kmer_pos_t);

// Use this to keep structs;
typedef struct {
  int x, y, nskip;
} segment_t;


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


/// utility function;
SEXP mk_strsxp(const char **words, size_t n){
  SEXP words_r = PROTECT(allocVector(STRSXP, n));
  for(size_t i=0; i < n; ++i)
    SET_STRING_ELT(words_r, i, mkChar(words[i]));
  UNPROTECT(1);
  return(words_r);
}

// We need a finalise function to clear resources
// when the external pointer goes out of scope
static void finalise_khash_ptr(SEXP ptr_r){
  void *v_ptr = R_ExternalPtrAddr(ptr_r);
  if(!v_ptr)
    return;
  khash_ptr *ptr = (khash_ptr*)v_ptr;
  if(ptr->hash){
    clear_kmer_h(ptr->hash);
    kh_destroy(kmer_h, ptr->hash);
    ptr->hash = 0;
  }
}

size_t skip_n(const char *seq, size_t i){
  while(seq[i] && LC(seq[i]) == 'n')
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
  }
  return( i + j );
}

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

// returns the total number of kmers encountered
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
      if(opt_flag & (flags[OPT_POS] | flags[OPT_PAIRS]) == 0)
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


static const R_CallMethodDef callMethods[] = {
	      {"make_kmer_h_index", (DL_FUNC)&make_kmer_h_index, 3},
	      {"kmer_positions", (DL_FUNC)&kmer_positions, 2},
	      {NULL, NULL, 0}
};

void R_init_kmer_hash(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
