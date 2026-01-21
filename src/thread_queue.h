#ifndef _THREAD_QUEUE_H
#define _THREAD_QUEUE_H
#include <pthread.h>
#include <stdint.h>
// data structs useful for multithreading the counting of k-mers

// A circular buffer
// m:   the size allocated to the buffer
// end: new data elements are added after the last position
//       at buffer[ (end + 1) % m ] and end updated accordingly
// added: incremented whenever an element is added to the queue      
// beg; buffer[beg] is removed by a processing thread, and
//      beg <- (beg + 1) % m
// removed: incremented whenever an element is removed from the queue
// For an element to be added [added < (removed + m)] must be TRUE
// otherwise the feeding thread has to either increase the size of the
// the buffer or call pthread_condition_wait( cond )
//
// For an element to be removed, [removed < added] must be true.
//
// Break from the processing loop when 
//
// processing the kmer (i.e. putting it into a hash or and index of some sort
// is likely to take longer than obtaining it; for lots of reads and a k=20
// this seems to be something like 75 times longer. This means that we'd like to
// have at least 75 threads, but it means that the reading thread will need to
// wait for the completion of work by worker threads.
//
// We can either have a single queue for all threads, or one queue for each thread.
// The latter seems reasonable to me as it should result in less contention for resources.
//
// buffer: the buffer
// m : the size of the buffer
// n : the number of elements held in the queue, pending processing
//     should be incremented when elements are added and decremented when
//     elements are removed.
// i : the postion in the queue from which elements should be removed
//     elements should be added at buffer[(i + n) % m]
//     if n < m
// total: the total number of elements added to the queue
// done: if done && n == 0, then the loop should exit.
typedef struct {
  uint64_t *buffer;
  size_t m; // the size of the buffer
  size_t n;
  size_t i;
  size_t total;
  int *done;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
} kmer_queue;

void init_kmer_queue(kmer_queue *kq, size_t buf_size, int *done);
void free_kmer_queue(kmer_queue *kq);

void kmer_queue_add(kmer_queue *kq, uint64_t kmer);
//uint64_t kmer_queue_remove(int *err);

#endif
