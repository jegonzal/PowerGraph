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


#ifndef GRAPHLAB_SYNCHRONIZED_QUEUE_HPP
#define GRAPHLAB_SYNCHRONIZED_QUEUE_HPP


#include <queue>

#include <graphlab/parallel/pthread_tools.hpp>
namespace graphlab {
  /// \ingroup util_internal
  template <typename T>
  class synchronized_queue {

  public:
    synchronized_queue() { };
    ~synchronized_queue() { };

    void push(const T &item) {
      _queuelock.lock();
      _queue.push(item);
      _queuelock.unlock();
    }
  
    bool safepop(T * ret) {
      _queuelock.lock();
      if (_queue.size() == 0) {
        _queuelock.unlock();

        return false;
      }
      *ret = _queue.front();
      _queue.pop();
      _queuelock.unlock();
      return true;
    }
  
    T pop() {
      _queuelock.lock();
      T t = _queue.front();
      _queue.pop();
      _queuelock.unlock();
      return t;
    }
  
    size_t size() const{
      return _queue.size();
    }
  private:
    std::queue<T> _queue;
    spinlock _queuelock;
  };

}
#endif

