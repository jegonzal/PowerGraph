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

#ifndef GRAPHLAB_SYNCHRONIZED_MULTIQUEUE_HPP
#define GRAPHLAB_SYNCHRONIZED_MULTIQUEUE_HPP
#include <queue>

#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/random.hpp>


namespace graphlab {
  /// \ingroup util_internal
  template <typename T>
  class synchronized_multiqueue {

  public:
    synchronized_multiqueue(size_t numqueues_ = 8) {
      numqueues = numqueues_;
      queue.resize(numqueues);
      locks.resize(numqueues);
      num_queues_full.value = 0;
    };
    ~synchronized_multiqueue() { };

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
      //  M.D. Mitzenmacher The Power of Two Choices in Randomized
      // Load Balancing (1991)
      // http://www.eecs.harvard.edu/~michaelm/postscripts/mythesis.pdf

      // pushing without a hint. push into the smaller of 2 selections
      if (numqueues == 1) { 
        push(item, 0);
        return;
      }
      size_t r1 = random::uniform<size_t>(0, numqueues-1);
      size_t r2 = random::uniform<size_t>(0, numqueues-2);
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


    // static inline size_t rand_int(size_t q) {
    //   //     double v = thread::rand01();
    //   //     size_t r = *reinterpret_cast<size_t*>(&v); 
    //   //     return (r >> 8) % q;
    //   return  size_t(std::floor(thread::rand01() * q));
    // }

    std::vector<std::queue<T> > queue;
    std::vector<spinlock> locks;
    size_t numqueues;
    atomic<size_t> num_queues_full;
  };

}
#endif

