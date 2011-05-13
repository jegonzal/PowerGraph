/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
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
        the blocking queue is essentially destroyed and can no longer be used.
    */
    inline void stop_blocking() {
      m_mutex.lock();
      m_alive = false;
      m_conditional.broadcast();
      m_empty_conditional.broadcast();
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
     * Returns true on success
     * Returns false if the queue is no longer alove
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

    /**
     * Causes any threads currently blocking on a dequeue to wake up
     */
    void signal() {
      m_mutex.lock();
      m_conditional.broadcast();
      m_mutex.unlock();
    }
    void signal_blocking_empty() {
      m_mutex.lock();
      m_empty_conditional.broadcast();
      m_mutex.unlock();
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

