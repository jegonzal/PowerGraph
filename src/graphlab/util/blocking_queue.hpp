/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef GRAPHLAB_BLOCKING_QUEUE_HPP
#define GRAPHLAB_BLOCKING_QUEUE_HPP



#include <list>
#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** 
    * \ingroup util
    * \brief Implements a blocking queue useful for producer/consumer models
    */
  template<typename T>
  class blocking_queue {
  private:
    
    typedef typename std::list<T> queue_type;

    bool m_alive;
    queue_type m_queue;
    mutex m_mutex;
    conditional m_conditional;
    conditional m_empty_conditional;
    
    /**
     * Causes any threads currently blocking on a dequeue to wake up
     * and evaluate the state of the queue. If the queue is empty,
     * the threads will return back to sleep immediately. If the queue
     * is destroyed through stop_blocking, all threads will return. 
     */
    void signal() {
      m_mutex.lock();
      m_conditional.broadcast();
      m_mutex.unlock();
    }
    
    

    /**
     * Causes any threads blocking on "wait_until_empty()" to wake
     * up and evaluate the state of the queue. If the queue is not empty,
     * the threads will return back to sleep immediately. If the queue
     * is empty, all threads will return.
    */
    void signal_blocking_empty() {
      m_mutex.lock();
      m_empty_conditional.broadcast();
      m_mutex.unlock();
    }    


  public:
    
    //! creates a blocking queue
    blocking_queue() : m_alive(true) { }
    
    //! Add an element to the blocking queue
    inline void enqueue(const T& elem) {
      m_mutex.lock();
      m_queue.push_back(elem);
      // Signal threads waiting on the queue
      m_conditional.signal();
      m_mutex.unlock();
    }

    /**
     * Blocks until an element is available in the queue 
     * or until stop_blocking() is called.
     * The return value is a pair of <T value, bool success>
     * If "success" if set, then "value" is valid and 
     * is an element popped from the queue.
     * If "success" is false, stop_blocking() was called 
     * and the queue has been destroyed.
     */
    inline std::pair<T, bool> dequeue() {

      m_mutex.lock();
      T elem = T();
      bool success = false;
      // Wait while the queue is empty and this queue is alive
      while(m_queue.empty() && m_alive) {
        m_conditional.wait(m_mutex);
      }
      // An element has been added or a signal was raised
      if(!m_queue.empty()) {
        success = true;
        elem = m_queue.front();
        m_queue.pop_front();
        if (m_queue.empty()) {
          m_empty_conditional.signal();
        }
      } 
      m_mutex.unlock();

      return std::make_pair(elem, success);
    }

    /**
    * Returns an element if the queue has an entry.
    * returns [item, false] otherwise.
    */
    inline std::pair<T, bool> try_dequeue() {
      m_mutex.lock();
      T elem = T();
      // Wait while the queue is empty and this queue is alive
      if (m_queue.empty() || m_alive == false) {
        m_mutex.unlock();
        return std::make_pair(elem, false);
      }
      else {
        elem = m_queue.front();
        m_queue.pop_front();
        if (m_queue.empty()) {
          m_empty_conditional.signal();
        }
      }
      m_mutex.unlock();

      return std::make_pair(elem, true);
    }

    //! Returns true if the queue is empty
    inline bool empty() { 
      m_mutex.lock();
      bool res = m_queue.empty();
      m_mutex.unlock();
      return res;
    }

    /** Wakes up all threads waiting on the queue whether 
        or not an element is available. Once this function is called,
        all existing and future dequeue operations will return with failure.
        Note that there could be elements remaining in the queue after 
        stop_blocking() is called. 
    */
    inline void stop_blocking() {
      m_mutex.lock();
      m_alive = false;
      m_conditional.broadcast();
      m_empty_conditional.broadcast();
      m_mutex.unlock();
    }

    /**
      Resumes operation of the blocking_queue. Future calls to
      dequeue will proceed as normal.
    */
    inline void start_blocking() {
      m_mutex.lock();
      m_alive = true;
      m_mutex.unlock();
    }
    
    //! get the current size of the queue
    inline size_t size() {
      m_mutex.lock();
      size_t size = m_queue.size();
      m_mutex.unlock();
      return size;
    }

    /**
     * The conceptual "reverse" of dequeue().
     * This function will block until the queue becomes empty, or 
     * until stop_blocking() is called.
     * Returns true on success. 
     * Returns false if the queue is no longer alive
    */
    bool wait_until_empty() {
      m_mutex.lock();
      // if the queue still has elements in it while I am still alive, wait
      while (m_queue.empty() == false && m_alive == true) {
        m_empty_conditional.wait(m_mutex);
      }
      m_mutex.unlock();
      // if I am alive, the queue must be empty. i.e. success
      // otherwise I am dead
      return m_alive;
    }
    
    ~blocking_queue() {
      m_alive = false;
      signal();
      signal_blocking_empty();
    }    
  }; // end of blocking_queue class
  

} // end of namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif

