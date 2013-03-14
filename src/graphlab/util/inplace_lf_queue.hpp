#ifndef GRAPHLAB_INPLACE_LOCKFREE_QUEUE_HPP
#define GRAPHLAB_INPLACE_LOCKFREE_QUEUE_HPP
#include <stdint.h>
#include <cstring>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/atomic_ops.hpp>
#include <utility>
namespace graphlab {

/*
 * A lock free queue where each element is a byte sequence,
 * where the first 8 bytes can be used for a next pointer.
 *
 * head is the head of the queue. Always sentinel.
 * tail is current last element of the queue.
 * completed is the last element that is completely inserted.
 * There can only be one thread dequeueing.
 *
 * On dequeue_all, the dequeu-er should use get_next() to get the
 * next element in the list. If get_next() returns NULL, it should spin
 * until not null, and quit only when end_of_dequeue_list() evaluates to true
 */
class inplace_lf_queue {
 public:
   inline inplace_lf_queue():head(sentinel),tail(sentinel) {
     for (size_t i = 0;i < sizeof(size_t); ++i) sentinel[i] = 0;
   }

   void enqueue(char* c);

   char* dequeue_all();

   static inline char* get_next(char* ptr) {
     return *(reinterpret_cast<char**>(ptr));
   }

   static inline char** get_next_ptr(char* ptr) {
     return reinterpret_cast<char**>(ptr);
   }

   inline const bool end_of_dequeue_list(char* ptr) {
     return ptr == sentinel;
   }

 private:

   char sentinel[sizeof(size_t)];
   char* head;
   char* tail;
};


} // namespace graphlab

#endif
