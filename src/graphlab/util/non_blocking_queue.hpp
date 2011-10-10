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


#ifndef GRAPHLAB_NON_BLOCKING_QUEUE_HPP
#define GRAPHLAB_NON_BLOCKING_QUEUE_HPP



#include <list>
#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** 
    * \ingroup util
    * \brief Implements a blocking queue useful for producer/consumer models
    */
  template<typename T>
  class non_blocking_queue {
  
   private:
    
    typedef typename std::list<T> queue_type;

    queue_type m_queue;
    simple_spinlock m_mutex;
    
   public:
    
    //! creates a blocking queue
    non_blocking_queue() { }
    
    //! Add an element to the blocking queue
    inline void enqueue(const T& elem) {
      m_mutex.lock();
      m_queue.push_back(elem);
      m_mutex.unlock();
    }


    void begin_critical_section() {
      m_mutex.lock();
    }

    void end_critical_section() {
     m_mutex.unlock();
    }

    inline std::pair<T, bool> get_front() {
      m_mutex.lock();
      T elem = T();
      // Wait while the queue is empty and this queue is alive
      if (m_queue.empty()) {
        m_mutex.unlock();
        return std::make_pair(elem, false);
      }
      else {
        elem = m_queue.front();
      }
      m_mutex.unlock();
      return std::make_pair(elem, true);
    }


    inline std::pair<T, bool> try_dequeue() {
      m_mutex.lock();
      T elem = T();
      // Wait while the queue is empty and this queue is alive
      if (m_queue.empty()) {
        m_mutex.unlock();
        return std::make_pair(elem, false);
      }
      else {
        elem = m_queue.front();
        m_queue.pop_front();
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
   
    ~non_blocking_queue() {
    }    
  }; // end of blocking_queue class
  

} // end of namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif

