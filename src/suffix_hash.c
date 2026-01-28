#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "suffix_hash.h"
#include "thread_queue.h"

suffix_hash init_suffix_hash(uint32_t prefix_bits, uint32_t suffix_bits, size_t max_size){
  suffix_hash kt;
  uint32_t total_bits = prefix_bits + suffix_bits;
  kt.suffix_bits = suffix_bits <= SH_MAX_SUF_BITS ? suffix_bits : SH_MAX_SUF_BITS;
  prefix_bits = total_bits - suffix_bits;
  kt.kmer_mask = (total_bits >= 64) ? ~(size_t)0 : ((size_t)1 << total_bits) - 1;
  kt.suffix_mask = ((size_t)1 << suffix_bits) - 1;
  kt.prefix_bits = prefix_bits;
  kt.prefixes = 0;
  kt.prefix_n = ((size_t)1 << prefix_bits);
  kt.max_size = max_size;
  kt.allocated = 0;
  kt.max_count = 0;
  kt.max_count_kmer = 0;
  if(sizeof(void*) * kt.prefix_n <= max_size)
    kt.prefixes = calloc(kt.prefix_n, sizeof(void*));
  return(kt);
}

void init_suffix_hash_p(suffix_hash *sh, uint32_t prefix_bits, uint32_t suffix_bits, size_t max_size){
  uint32_t total_bits = prefix_bits + suffix_bits;
  sh->suffix_bits = suffix_bits <= SH_MAX_SUF_BITS ? suffix_bits : SH_MAX_SUF_BITS;
  prefix_bits = total_bits - suffix_bits;
  sh->kmer_mask = (total_bits >= 64) ? ~(size_t)0 : ((size_t)1 << total_bits) - 1;
  sh->suffix_mask = ((size_t)1 << suffix_bits) - 1;
  sh->prefix_bits = prefix_bits;
  sh->prefixes = 0;
  sh->prefix_n = ((size_t)1 << prefix_bits);
  sh->max_size = max_size;
  sh->allocated = 0;
  sh->max_count = 0;
  sh->max_count_kmer = 0;
  if(sizeof(void*) * sh->prefix_n <= max_size)
    sh->prefixes = calloc(sh->prefix_n, sizeof(void*));
}


void free_suffix_hash(suffix_hash *kt){
  if(!kt->prefixes)
    return;
  for(size_t i=0; i < kt->prefix_n; ++i){
    if(kt->prefixes[i])
      kh_destroy(kcount, kt->prefixes[i]);
  }
  free(kt->prefixes);
}

// returns negative values for errors
// else, the actual count of the kmer
int sh_add_kmer(suffix_hash *sh, uint64_t kmer){
  kmer &= sh->kmer_mask;
  size_t prefix_i = kmer >> sh->suffix_bits;
  size_t suffix = kmer & sh->suffix_mask;
  if(prefix_i >= sh->prefix_n)
    return(-1);
  if(!sh->prefixes[prefix_i]){
    sh->prefixes[prefix_i] = kh_init(kcount);
    if(!sh->prefixes[prefix_i]){
      printf("Failed to initialise new hash for prefix: %ld suffix: %ld   kmer: %ld\n", prefix_i, suffix, kmer);
      printf("previously allocated: %ld\n", sh->allocated);
      printf("prefix / suffix bits %d / %d\n", sh->prefix_bits, sh->suffix_bits);
      return(-2);
    }
    sh->allocated++;
  }
  khash_t(kcount) *hash = sh->prefixes[prefix_i];
  int ret;
  khiter_t k = kh_get(kcount, hash, (uint32_t)suffix);
  if(k == kh_end(hash) || !kh_exist(hash, k)){
    k = kh_put(kcount, hash, suffix, &ret);
    kh_value(hash, k) = 1;
  }else{
    kh_value(hash, k)++;
  }
  uint32_t count = kh_value(hash, k);
  if(count > sh->max_count){
    sh->max_count = count;
    sh->max_count_kmer = kmer;
  }
  return((int)count);
}

uint32_t sh_kmer_count(suffix_hash *sh, uint64_t kmer){
  kmer &= sh->kmer_mask;
  size_t p_i = kmer >> sh->suffix_bits;
  size_t s_i = kmer & sh->suffix_mask;
  if(p_i >= sh->prefix_n || !sh->prefixes[p_i])
    return(0);
  khash_t(kcount) *hash = sh->prefixes[p_i];
  khiter_t k = kh_get(kcount, hash, (uint32_t)s_i);
  if(k == kh_end(hash) || !kh_exist(hash, k))
    return(0);
  return(kh_value(hash, k));
}

size_t sh_count_spectrum(suffix_hash *sh, double *counts, uint32_t counts_n){
  uint32_t max_count = counts_n - 1;
  size_t n = 0;
  for(size_t i=0; i < sh->prefix_n; ++i){
    if(!sh->prefixes[i])
      continue;
    khash_t(kcount) *hash = sh->prefixes[i];
    khiter_t k;
    for(k=kh_begin(hash); k != kh_end(hash); ++k){
      if(!kh_exist(hash, k)) continue;
      uint32_t count = kh_value(hash, k);
      counts[ count >= max_count ? max_count : count ]++;
      if(count > 0)
	++n;
    }
  }
  return(n);
}

int init_suffix_hash_n(suffix_hash_n* sh, uint32_t counts_n,
			uint32_t k, uint32_t prefix_bits, uint32_t suffix_bits){
  uint32_t total_bits = prefix_bits + suffix_bits;
  if(total_bits != k * 2)
    return(-1);
  if(total_bits > SH_KMAX_BITS)
    return(-2);
  if(prefix_bits > SH_MAX_PRE_BITS)
    return(-3);
  if(suffix_bits > SH_MAX_SUF_BITS)
    return(-4);
  if(counts_n > 4)
    return(-5);
  sh->k = k;
  sh->kmer_mask = (1ULL << total_bits) - 1;
  sh->suffix_mask = (1ULL << sh->suffix_bits);
  sh->prefix_bits = prefix_bits;
  sh->suffix_bits = suffix_bits;
  sh->counts_n = counts_n;
  sh->prefix_n = 1ULL << prefix_bits;
  sh->prefixes = calloc(sh->prefix_n, sizeof(void**));
  sh->kmer_counts = calloc(counts_n, sizeof(uint64_t));
  return(1);
}

void free_suffix_hash_n(suffix_hash_n *sh){
  for(size_t i=0; i < sh->prefix_n; ++i){
    switch(sh->counts_n){
    case 1: // khash_t(kcount)
      kh_destroy(kcount, (khash_t(kcount)*)sh->prefixes[i]);
      break;
    case 2:
      kh_destroy(kcount_2, (khash_t(kcount_2)*)sh->prefixes[i]);
      break;
    case 3:
      kh_destroy(kcount_3, (khash_t(kcount_3)*)sh->prefixes[i]);
      break;
    case 4:
      kh_destroy(kcount_4, (khash_t(kcount_4)*)sh->prefixes[i]);
      break;
    default:
      printf("bugger; unimplemented suffix_hash_n destruction not possible\n");
      return;
    }
  }
  free(sh->prefixes);
  free(sh->kmer_counts);
}

int sh_n_add_kmer(suffix_hash_n *sh, uint32_t source, uint64_t kmer){
  if(source >= sh->counts_n)
    return(-1);
  kmer &= sh->kmer_mask;
  uint64_t prefix_i = kmer >> sh->suffix_bits;
  uint32_t suffix = (uint32_t)(kmer & sh->suffix_mask);
  if(prefix_i >= sh->prefix_n)
    return(-2);

  void *hash_p = sh->prefixes[prefix_i];
  if(!sh->prefixes[prefix_i]){
    switch(sh->counts_n){
    case 1:
      hash_p = kh_init(kcount);
      break;
    case 2:
      hash_p = kh_init(kcount_2);
      break;
    case 3:
      hash_p = kh_init(kcount_3);
      break;
    case 4:
      hash_p = kh_init(kcount_4);
      break;
    default:
      printf("sh_n_add_kmer; kcount_n too high\n");
      return(-3);
    }
  }
    // This is rather ugly, but can't come up with  better way for now.
  khiter_t k;
  int ret;
  switch(sh->counts_n){
  case 1:
    khash_t(kcount) *hash = (khash_t(kcount)*)hash_p;
    k = kh_get(kcount, hash, suffix);
    if(k == kh_end(hash) || !kh_exist(hash, k)){
      k = kh_put(kcount, hash, suffix, &ret);
      kh_value(hash, k) = 1;
    }else{
      kh_value(hash, k)++;
    }
    break;
  case 2:
    khash_t(kcount_2) *hash_2 = (khash_t(kcount_2)*)hash_p;
    k = kh_get(kcount_2, hash_2, suffix);
    if(k == kh_end(hash_2) || !kh_exist(hash_2, k)){
      k = kh_put(kcount_2, hash_2, suffix, &ret);
      memset(&kh_value(hash_2, k), 0, sizeof(uint32_t) * sh->counts_n);
      kh_value(hash_2, k).n[source] = 1;
    }else{
      kh_value(hash_2, k).n[source]++;
    }
    break;
  case 3:
    khash_t(kcount_3) *hash_3 = (khash_t(kcount_3)*)hash_p;
    k = kh_get(kcount_3, hash_3, suffix);
    if(k == kh_end(hash_3) || !kh_exist(hash_3, k)){
      k = kh_put(kcount_3, hash_3, suffix, &ret);
      memset(&kh_value(hash_3, k), 0, sizeof(uint32_t) * sh->counts_n);
      kh_value(hash_3, k).n[source] = 1;
    }else{
      kh_value(hash_3, k).n[source]++;
    }
    break;
  case 4:
    khash_t(kcount_4) *hash_4 = (khash_t(kcount_4)*)hash_p;
    k = kh_get(kcount_4, hash_4, suffix);
    if(k == kh_end(hash_4) || !kh_exist(hash_4, k)){
      k = kh_put(kcount_4, hash_4, suffix, &ret);
      memset(&kh_value(hash_4, k), 0, sizeof(uint32_t) * sh->counts_n);
      kh_value(hash_4, k).n[source] = 1;
    }else{
      kh_value(hash_4, k).n[source]++;
    }
    break;
  default:
    printf("sh_n_add_kmer; kcount_n too high\n");
    return(-4);
  }
  return(1);
}

size_t sh_count_spectrum_nc(suffix_hash_n *sh, uint32_t *counts, uint32_t counts_l,
			    uint32_t comb_in, uint32_t inner){
  return(0);
}


/* void *sh_intersect_thread( void* args ){ */
/*   return(args); */
/* } */

/* suffix_hash_2 sh_intersect(suffix_hash *sh1, suffix_hash *sh2, uint32_t thread_n){ */
/*   suffix_hash_2 sh_pair; */
/*   // init and do something */
/*   return(sh_pair); */
/* } */


// NOTE: This approach of multithreading did not work well;
//       in fact it only made things slower. I should probably
//       remove it; however, it may be that modifying in some
//       way or other might work OK.
// this is a bit wasteful. We could do with a simple pointer as the
// only thing that is different between threads is the kmer_queue
void init_thread_args(suffix_hash_t_args *args,
		      uint32_t prefix_bits, uint32_t suffix_bits,
		      khash_t(kcount) **hashes, kmer_queue *queue,
		      int *done, size_t thread_i){
  args->prefix_bits = prefix_bits;
  args->suffix_bits = suffix_bits;
  args->suffix_mask = (1ULL << suffix_bits) - 1;
  //  args->prefix_n = (1ULL << prefix_bits); // maybe not needed?
  args->hashes = hashes;
  args->queue = queue;
  args->done = done;
  args->thread_i = thread_i;
  args->hashes_allocated = 0;
  args->kmer_n = 0;
}

// the mutex and condition in the args are expected to already be initialised.
// WARNING:
//    this function assumes that the thread will have ownership of the
void *sh_process_kmers(void *thread_args){
  suffix_hash_t_args *args = (suffix_hash_t_args*)thread_args;
  // lock the mutex, and call condition wait
  pthread_mutex_lock( &args->queue->mutex );
  while(1){
    pthread_cond_wait( &args->queue->cond, &args->queue->mutex );
    // I have data to process; at this point I hold the lock;
    while(args->queue->n){
      //      pthread_mutex_lock( &args->queue->mutex );
      //printf("%ld\t%ld %ld %ld done: %d\n", args->thread_i, args->queue->i, args->queue->n, args->queue->m, *args->done);
      uint64_t kmer = args->queue->buffer[ args->queue->i ];
      args->queue->n--;
      args->queue->i++;
      args->queue->i %= args->queue->m;
      //      pthread_mutex_unlock(&args->queue->mutex);
      // It is possible that the reader thread is waiting to add a k-mer
      // I'm not sure whether it is better to signal immediately, or to signal
      // after the k-mer has been processed. Signalling here means the reader can
      // get back to reading faster, but may lead to a lot of un-necessary signalling.
      pthread_cond_signal(&args->queue->cond);
      size_t prefix_i = kmer >> args->suffix_bits;
      uint32_t suffix = kmer & args->suffix_mask;
      if(!args->hashes[prefix_i]){
	args->hashes[prefix_i] = kh_init(kcount);
	args->hashes_allocated++;
      }
      khash_t(kcount) *hash = args->hashes[prefix_i];
      khiter_t k = kh_get(kcount, hash, suffix);
      int ret = 0;
      if(k == kh_end(hash)){
	k = kh_put(kcount, hash, suffix, &ret);
	kh_value(hash, k) = 1;
	args->kmer_n++;
      }else{
	kh_value(hash, k)++;
      }
      pthread_mutex_lock( &args->queue->mutex );
    }
    if(*args->done && args->queue->n == 0){
      pthread_mutex_unlock( &args->queue->mutex );
      break;
    }
  }
  return(thread_args);
}

void init_suffix_hash_mt(suffix_hash_mt *sh,
			 uint32_t prefix_bits, uint32_t suffix_bits,
			 size_t nthreads, size_t queue_buffer_size){
  uint32_t total_bits = prefix_bits + suffix_bits;
  total_bits = total_bits > SH_KMAX_BITS ? SH_KMAX_BITS : total_bits;
  sh->suffix_bits = suffix_bits <= SH_MAX_SUF_BITS ? suffix_bits : SH_MAX_SUF_BITS;
  sh->prefix_bits = total_bits - suffix_bits;
  sh->prefix_bits = prefix_bits;
  sh->kmer_mask = (1ULL << (prefix_bits + suffix_bits));
  sh->prefix_n = (1ULL << prefix_bits);
  sh->done = 0;
  sh->hashes = calloc(sh->prefix_n, sizeof(khash_t(kcount)*));
  sh->nthreads = nthreads;
  sh->t_args = calloc(nthreads, sizeof(suffix_hash_t_args));
  sh->queues = calloc(nthreads, sizeof(kmer_queue));
  sh->threads = calloc(nthreads, sizeof(pthread_t));
  for(size_t i=0; i < nthreads; ++i){
    init_kmer_queue(sh->queues + i, queue_buffer_size, &sh->done);
    init_thread_args(sh->t_args + i, prefix_bits, suffix_bits,
				    sh->hashes, &sh->queues[i], &sh->done, i);
    pthread_create(&sh->threads[i], NULL, &sh_process_kmers,
		   &sh->t_args[i]);
  }
}

// This assumes that all threads have already been joined.
// if that hasn't been done then I'm not sure what will happen.
void free_suffix_hash_mt(suffix_hash_mt *sh){
  for(size_t i=0; i < sh->prefix_n; ++i){
    if(sh->hashes[i])
      kh_destroy(kcount, sh->hashes[i]);
  }
  // it seems that there is not any need to destroy pthreads?
  free(sh->hashes);
  free(sh->threads);
  free(sh->queues);
}

// return 0 on success
int suffix_hash_mt_add(suffix_hash_mt* sh, uint64_t kmer){
  size_t prefix = (kmer >> sh->suffix_bits);
  if(prefix >= sh->prefix_n)
    return(-1);
  size_t thread_i = prefix % sh->nthreads;
  kmer_queue_add(&sh->queues[thread_i], kmer);
  return(0);
}

/* void suffix_hash_mt_join(suffix_hash_mt *sh){ */
/*   sh->done = 1; */
/*   for(size_t i=0; i < sh->nthreads; ++i){ */
/*     pthread_cond_signal(&sh->queues[i].cond); */
/*     pthread_join(sh->threads[i], NULL); */
/*     printf("thread %ld joined. Allocated %ld. kmers: %ld\n", i, sh->t_args[i].hashes_allocated, sh->t_args[i].kmer_n); */
/*   } */
/* } */
