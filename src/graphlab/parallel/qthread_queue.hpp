#ifndef QTHREAD_TASK_QUEUE_HPP
#define QTHREAD_TASK_QUEUE_HPP
#include <deque>
#include <graphlab/parallel/qthread_tools.hpp>
namespace graphlab {

/**
 * This implements a nonblocking queue
 * with qthread routines. Note that this class 
 * is only thread safe when used inside qthread threads.
 * This class cannot be used safely using regular threads
 * and may fault in interesting ways.
 */
template <typename T>
class qthread_queue {
 private:
  typedef typename std::deque<T> queue_type;
  bool m_alive;
  qthread_mutex m_lock;
  queue_type m_queue;
   
 public:
  //! creates a blocking queue
  qthread_queue() : m_alive(true) { }

  void enqueue(const T& elem) {
    m_mutex.lock();
    m_queue.push_back(elem);
    // Signal threads waiting on the queue
    m_mutex.unlock();
  }


  //! Add an element to the blocking queue
  void enqueue_to_head(const T& elem) {
    m_mutex.lock();
    m_queue.push_front(elem);
    // Signal threads waiting on the queue
    m_mutex.unlock();
  }

  /**
   * Returns true if an element is returned and false
   * otherwise. The element is returned in *result if result != NULL
   */
  inline bool dequeue_nonblocking(T* result) {
    bool ret = false;
    m_mutex.lock();
    if(m_alive && !m_queue.empty()){ 
      if (result != NULL) (*result) = m_queue.front();
      m_queue.pop_front();
      ret = true;
    }
    m_mutex.unlock();
    return ret;
  }

  bool alive() const {
    return m_alive;
  }
  
  /** Wakes up all threads waiting on the queue whether 
    or not an element is available. Once this function is called,
    all existing and future dequeue operations will return with failure.
    Note that there could be elements remaining in the queue after 
    stop_blocking() is called. 
   */
  inline void stop_blocking() {
    m_alive = false;
  }
}; 

} // namespace graphlab
#endif
