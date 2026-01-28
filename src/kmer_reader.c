#include <stdio.h>
#include <zlib.h>
#include <float.h> // for DBL_MAX
#include <math.h>
#include "kmer_reader.h"

// -708 is approximately log(-DBL_MAX)
// it doesn't really matter so much, but it's a reasonable value to use.
// The other values here were obtained using R and copy and paste:
// log( 1 - 10^(-q/10) )
static const double q_to_ll[256] = {-708, -708, -708, -708, -708, -708, -708, -708,
				       -708, -708, -708, -708, -708, -708, -708, -708,
				       -708, -708, -708, -708, -708, -708, -708, -708,
				       -708, -708, -708, -708, -708, -708, -708, -708,
				       -708, -708,
				       -1.58147375340845, -0.996843044007847, -0.695524471332314, -0.507675873696745, -0.380130408066172, -0.289268187201608,
				       -0.222551515972833, -0.172556572913836, -0.134551960288774, -0.105360515657826, -0.0827653026691603, -0.0651741731994309,
				       -0.0514182741579622, -0.0406248442216418, -0.0321335740230858, -0.0254397275342231, -0.0201543647612323,
				       -0.0159758692466999, -0.0126691702086354, -0.0100503358535015, -0.00797499827851269, -0.00632956293111125,
				       -0.00502447389098515, -0.00398901726640658, -0.00316728822615736, -0.002515046511182, -0.00199725550255094,
				       -0.00158615046428019, -0.00125971852410645, -0.00100050033358353, -0.000794643880558492, -0.000631156481834621,
				       -0.000501312869928779, -0.000398186436251368, -0.000316277776560234, -0.000251220196302187, -0.000199546139503567,
				       -0.000158501880005453, -0.000125900466310563, -0.000100005000333347, -7.94359784262421e-05, -6.30977250676572e-05,
				       -5.01199793478725e-05, -3.98115095230173e-05, -3.16232766122775e-05, -2.51191797990654e-05, -1.99528222059213e-05,
				       -1.58490575202822e-05, -1.25893333632809e-05, -1.00000500002878e-05, -7.94331389524571e-06, -6.30959335027155e-06,
				       -5.01188489573212e-06, -3.98107963000364e-06, -3.16228266018417e-06, -2.51188958631412e-06, -1.99526430546292e-06,
				       -1.5848944484079e-06, -1.25892620425398e-06, -1.00000050002909e-06, -7.94328550222186e-07, -6.30957543536599e-07,
				       -5.01187359197959e-07, -3.98107249807448e-07, -3.16227815972964e-07, -2.51188674655675e-07, -1.99526251442212e-07,
				       -1.58489331783564e-07, -1.25892549094037e-07, -1.00000004947365e-07, -7.94328265847126e-08, -6.30957364055222e-08,
				       -5.01187246496094e-08, -3.98107178487233e-08, -3.16227771194999e-08, -2.51188646373612e-08, -1.99526233083297e-08,
				       -1.58489320702119e-08, -1.2589254162895e-08, -1.00000001002476e-08, -7.9432823967449e-09, -6.30957343919953e-09,
				       -5.01187237413051e-09, -3.98107170244991e-09, -3.16227766694999e-09, -2.51188647975196e-09, -1.99526229071369e-09,
				       -1.58489321792216e-09, -1.25892540915748e-09, -9.99999972218069e-10, -7.94328270141873e-10, -6.30957397639622e-10,
				       -5.01187202976012e-10, -3.9810721394071e-10, -3.16227710733847e-10, -2.51188625486804e-10, -1.99526284403372e-10,
				       -1.58489332781411e-10, -1.25892518639967e-10, -1.00000008279037e-10, -7.94327936791035e-11, -6.30957508482776e-11,
				       -5.01186869796069e-11, -3.98107102847087e-11, -3.16228154778057e-11, -2.51189069547621e-11, -1.9952595131855e-11,
				       -1.58488777658596e-11, -1.25892629655137e-11, -1.00000008274537e-11, -7.94331267431716e-12, -6.3096194935696e-12,
				       -5.0118798000779e-12, -3.98103772170881e-12, -3.16224824104483e-12, -2.51187959321757e-12, -1.99529281985832e-12,
				       -1.58484336765367e-12, -1.25888188762326e-12, -9.99977878280379e-13, -7.94364574119615e-13, -6.30939744894675e-13,
				       -5.01154673315921e-13, -3.9812597663066e-13, -3.16191517413295e-13, -2.51243470472704e-13, -1.99507077525161e-13,
				       -1.58539847916485e-13, -1.25899290992501e-13, -1.00031094518732e-13, -7.93809462607018e-14, -6.30606677987109e-14,
				       -5.00710584105958e-14, -3.98570065840439e-14, -3.16413562018175e-14, -2.50910403565289e-14, -1.9984014443253e-14,
				       -1.58761892521399e-14, -1.25455201782643e-14, -9.99200722162646e-15, -7.99360577730116e-15, -6.32827124036341e-15,
				       -4.99600361081322e-15, -3.99680288865057e-15, -3.10862446895044e-15, -2.55351295663786e-15, -1.99840144432528e-15,
				       -1.55431223447522e-15, -1.22124532708767e-15, -9.99200722162641e-16, -7.7715611723761e-16, -6.66133814775094e-16,
				       -5.55111512312578e-16, -4.44089209850063e-16, -3.33066907387547e-16, -2.22044604925031e-16, -2.22044604925031e-16,
				       -1.11022302462516e-16, -1.11022302462516e-16, -1.11022302462516e-16, -1.11022302462516e-16, -1.11022302462516e-16,
				       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


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

int kr_add_kmer(kmer_reader *kr, uint64_t kmer_f, uint64_t kmer_r,
		size_t *word_count, size_t *kmer_count){
  uint64_t kmer = kmer_f < kmer_r ? kmer_f : kmer_r;
  size_t prefix_i = kmer >> kr->suffix_bits;
  int ins_ret = 0;
  if(!kr->restricted || (prefix_i % kr->thread_n) == kr->thread_i){
    ins_ret = sh_add_kmer(kr->sh, kmer);
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
  //uint64_t off_fwd = 0;
  //uint64_t off_rev = 0;
  //uint64_t offset = 0;
  //  uint64_t offset_rc = 0;
  //  uint64_t kmer = 0;
  //  uint64_t mask = (1ULL << (kr->k*2)) - 1;
  size_t word_count = 0;
  size_t kmer_count = 0;
  unsigned char *qual;
  unsigned char *seq;
  //  uint32_t rc_shift = 64 - kr->k * 2;
  //   int ins_ret = 0;
  //  uint64_t prefix_i = 0;
  kmer_iterator km_it;
  kmer_iterator_init(&km_it, kr->k, kr->min_q);
  while( (l = kseq_read(ks)) >= 0 && read_n < kr->max_reads ){
    ++read_n;
    if(l <= kr->k)
      continue;
    // The following would be better off as a separate function:
    // seq_to_kmer( struct seq_k )
    // where seq_k, contains the kseq pointer as well as the state (i.e. position
    // and if initialised)
    //    size_t i = 0;
    qual = (unsigned char*)(ks->qual.l == ks->seq.l ? ks->qual.s : 0);
    seq = (unsigned char*)ks->seq.s;
    uint64_t kmer_f=0;
    uint64_t kmer_r=0;
    // using these iterators makes the code quite a bit slower;
    // I should consider to write a function that procsses one read
    // using the ugly code below... but that's for later.
    if(!kmer_iterator_begin(&km_it, seq, qual, &kmer_f, &kmer_r))
      continue;
    if( kr_add_kmer(kr, kmer_f, kmer_r, &word_count, &kmer_count) < 0 )
      continue;
    while(kmer_iterator_next(&km_it, &kmer_f, &kmer_r)){
      if( kr_add_kmer(kr, kmer_f, kmer_r, &word_count, &kmer_count) < 0 )
	break;
    }
  }
  /*   while(seq[i]){ */
  /*   // pass any potential Ns and low quality regions */
  /*     i = init_kmer_qual_2(seq, qual, kr->min_q, i, &offset, &offset_rc, kr->k); */
  /*     if(!seq[i]) */
  /* 	break; */
  /*     off_fwd = offset & mask; */
  /*     off_rev = (offset_rc >> rc_shift) & mask; */
  /*     kmer = off_fwd < off_rev ? off_fwd : off_rev; */
  /*     prefix_i = kmer >> kr->suffix_bits; */
  /*     if(!kr->restricted || (prefix_i % kr->thread_n) == kr->thread_i){ */
  /* 	ins_ret = sh_add_kmer(kr->sh, kmer); */
  /* 	if(ins_ret < 0) */
  /* 	  break; */
  /* 	++word_count; */
  /* 	if(ins_ret == 1) */
  /* 	  ++kmer_count; */
  /*     } */
  /*     while(seq[i] && LC(seq[i]) != 'n' && (!qual || qual[i] > kr->min_q)){ */
  /* 	offset = UPDATE_OFFSET(offset, seq[i]); */
  /* 	offset_rc = UPDATE_OFFSET_RC(offset_rc, seq[i]); */
  /* 	++i; */
  /* 	off_fwd = mask & offset; */
  /* 	off_rev = mask & ( offset_rc >> rc_shift ); */
  /* 	kmer = off_fwd < off_rev ? off_fwd : off_rev; */
  /* 	prefix_i = kmer >> kr->suffix_bits; */
  /* 	if(!kr->restricted || (prefix_i % kr->thread_n) == kr->thread_i){ */
  /* 	  ins_ret = sh_add_kmer(kr->sh, kmer); */
  /* 	  if(ins_ret < 0) */
  /* 	    break; // this will only break out of the current region; not the whole sequence. */
  /* 	  ++word_count; */
  /* 	  if(ins_ret == 1) */
  /* 	    ++kmer_count; */
  /* 	} */
  /*     } */
  /*   } */
  /* } */
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

int seq_kmer_counts(const char* seq, size_t seq_l, int* counts, suffix_hash *sh, int k){
  uint64_t off_f = 0;
  uint64_t off_r = 0;
  uint64_t kmer_f = 0;
  uint64_t kmer_r = 0;
  uint64_t kmer = 0;
  uint32_t rc_shift = 64 - (k * 2);
  uint64_t mask = (1ULL << (k * 2)) - 1;
  memset(counts, 0, sizeof(int) * seq_l);
  if(k*2 != sh->suffix_bits + sh->prefix_bits)
    return(-1);
  size_t i = 0;
  while(seq[i]){
    if(i == 0 || LC(seq[i]) == 'n'){
      i = init_kmer_qual_2(seq, 0, 0, i, &off_f, &off_r, k);
      kmer_f = off_f & mask;
      kmer_r = off_r >> rc_shift;
      kmer = (kmer_f < kmer_r) ? kmer_f : kmer_r;
      counts[i-k] = sh_kmer_count(sh, kmer);
      if(!seq[i])
	break;
    }
    off_f = UPDATE_OFFSET(off_f, seq[i]);
    off_r = UPDATE_OFFSET_RC(off_r, seq[i]);
    kmer_f = off_f & mask;
    kmer_r = off_r >> rc_shift;
    kmer = (kmer_f < kmer_r) ? kmer_f : kmer_r;
    counts[i-k+1] = sh_kmer_count(sh, kmer);
    ++i;
  }
  return(1);
}
  
