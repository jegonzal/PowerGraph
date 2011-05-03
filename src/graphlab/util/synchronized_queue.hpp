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
