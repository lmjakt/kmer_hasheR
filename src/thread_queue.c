#include <stdlib.h>
#include <stdio.h>
#include "thread_queue.h"

void init_kmer_queue(kmer_queue *kq, size_t buf_size, int *done){
  kq->m = buf_size;
  kq->n = 0;
  kq->i = 0;
  kq->buffer = malloc(sizeof(uint64_t) * kq->m);
  kq->total = 0;
  kq->done = done;
  pthread_mutex_init( &kq->mutex, NULL );
  pthread_cond_init( &kq->cond, NULL );
}

void free_kmer_queue(kmer_queue *kq){
  kq->m = 0;
  free(kq->buffer); kq->buffer = 0;
  kq->n = kq->i = kq->m = 0;
  kq->total = 0;
  pthread_mutex_destroy(&kq->mutex);
  pthread_cond_destroy(&kq->cond);
}

// if 
void kmer_queue_add(kmer_queue *kq, uint64_t kmer){
  pthread_mutex_lock( &kq->mutex );
  while(1){
    if(kq->n < kq->m){
      kq->buffer[ (kq->n + kq->i) % kq->m ] = kmer;
      kq->n++;
      kq->total++;
      //      printf("added i: %ld  n: %ld  pos: %ld\n", kq->i, kq->n, (kq->n + kq->i) % kq->m);
      pthread_mutex_unlock(&kq->mutex);
      break;
    }
    pthread_cond_wait( &kq->cond, &kq->mutex );
  }
  pthread_cond_signal(&kq->cond);
}
