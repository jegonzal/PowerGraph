#include <graphlab/util/inplace_lf_queue.hpp>
namespace graphlab {
void inplace_lf_queue::enqueue(char* c) {
  // clear the next pointer
  (*get_next_ptr(c)) = NULL;
  // atomically,
  // swap(tail, c)
  // tail->next = c;
  char* prev = c;
  atomic_exchange(tail, prev);
  (*get_next_ptr(prev)) = c;
  asm volatile ("" : : : "memory");
}


char* inplace_lf_queue::dequeue_all() {
  // head is the sentinel
  char* ret_head = get_next(head);
  if (ret_head == NULL) return NULL;
  // now, the sentinel is not actually part of the queue.
  // by the time get_next(sentinel) is non-empty, enqueue must have completely
  // finished at least once, since the next ptr is only connected in line 11.
  // enqueue the sentinel. That will be the new head of the queue.
  // Anything before the sentinel is "returned". And anything after is part
  // of the queue
  enqueue(sentinel);

  // The last element in the returned queue
  // will point to the sentinel.
  return ret_head;
}



} // namespace graphlab
