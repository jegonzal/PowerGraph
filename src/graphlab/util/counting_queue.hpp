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


#ifndef GRAPHLAB_COUNTING_QUEUE_HPP
#define GRAPHLAB_COUNTING_QUEUE_HPP

#include <list>
#include <deque>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/blocking_queue.hpp>
#include <graphlab/logger/logger.hpp>


#include <graphlab/macros_def.hpp>

namespace graphlab {
   /** 
    * \ingroup util
    * \brief Implements a blocking queue with counter, useful for producer/consumer models.
    * Pop the head after N consumers have touched it.
    */

  template<typename T>
    class counting_queue :
      public blocking_queue<T> {
  public:
    typedef blocking_queue<T> base;
    using base::m_alive;
    using base::m_mutex;
    using base::m_conditional;
    using base::m_empty_conditional;
    using base::m_queue;
    using base::sleeping;
    using base::sleeping_on_empty;

  private:
      size_t limit;
      size_t num_access;

  public:
      counting_queue(size_t n) : 
        base(),
        limit(n),
        num_access(0) {}
      
      inline std::pair<T, bool> poll_till_pop() {
        m_mutex.lock();
        T elem = T();
        bool success =false;

        while(m_queue.empty() && m_alive) {
          sleeping++;
          m_conditional.wait(m_mutex);
          sleeping--;
        }
      
        if(!m_queue.empty()) {
          success = true;

          elem = m_queue.front();
          ++num_access;

          if (num_access == limit) {
            logger(LOG_DEBUG , "Counting queue pop.\n"); 
            m_queue.pop_front();
            num_access = 0;
            if (m_queue.empty() && sleeping_on_empty) {
              m_empty_conditional.signal();
            }
          }
        } 
        m_mutex.unlock();

        return std::make_pair(elem, success);
      }


  }; // end of counting_queue class
} // end namespace graphlab

#include <graphlab/macros_undef.hpp>
#endif
