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


#ifndef GRAPHLAB_PRIORITIZED_MULTIQUEUE_HPP
#define GRAPHLAB_PRIORITIZED_MULTIQUEUE_HPP
#include <queue>

#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

namespace graphlab {

  ///  \ingroup util_internal
  template <typename T>
  class prioritized_multiqueue {

  public:
    prioritized_multiqueue(size_t numqueues_ = 8) {
      numqueues = numqueues_;
      queue.resize(numqueues);
      locks.resize(numqueues);
      num_queues_full.value = 0;
    };
    ~prioritized_multiqueue() { };

    inline void push(const T &item, size_t queuehint) {
      locks[queuehint].lock();
      if (queue[queuehint].empty()) num_queues_full.inc();
      queue[queuehint].push(item);
      locks[queuehint].unlock();
    }
  
    inline bool safepop(T* ret, size_t queuehint) {
      locks[queuehint].lock();
      // if the queue has stuff in it
      if (queue[queuehint].size() == 0) {
        locks[queuehint].unlock();
        return false;
      }
      *ret = queue[queuehint].front();
      queue[queuehint].pop();
      if (queue[queuehint].empty()) num_queues_full.dec();
      locks[queuehint].unlock();
      return true;
    }
  
    /** Pushes an element into the queue without a hint. If you use this, you 
        really should lock. */ 
    void push(const T &item) {
      //  M.D. Mitzenmacher The Power of Two Choices in Randomized Load Balancing (1991)
      // http://www.eecs.harvard.edu/~michaelm/postscripts/mythesis.pdf

      // pushing without a hint. push into the smaller of 2 selections
      if (numqueues == 1) { 
        push(item, 0);
        return;
      }
      size_t r1 = rand_int(numqueues);
      size_t r2 = rand_int(numqueues-1);
      if (r2 >= r1) ++r2;
      size_t queuehint = (queue[r1].size() < queue[r2].size()) ? r1 : r2;
      push(item, queuehint);
    }
  
    size_t size() const{
      size_t count = 0;
      for (size_t i = 0; i < queue.size(); ++i) {
        count += queue[i].size();
      }
      return count;
    }

    bool empty() const {
      return (num_queues_full.value == 0);
    }
  private:
   
    static inline size_t rand_int(size_t q) {
      double v = random::rand01();
      size_t r = *reinterpret_cast<size_t*>(&v); 
      return (r >> 8) % q;
    }

    typedef mutable_queue<size_t, double> vertex_queue_type;

    std::vector<std::queue<T> > queue;
    std::vector<spinlock> locks;
    size_t numqueues;
    atomic<size_t> num_queues_full;
  };

}
#endif

