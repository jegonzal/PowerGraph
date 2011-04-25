#ifndef LAZY_DEQUE_HPP
#define LAZY_DEQUE_HPP
#include <deque>
#include <algorithm>
namespace graphlab {


/**
The lazy deque is a custom datastructure built with the following properties.
 - push_back is fast ( O(1) )
 - deletion from arbitrary locations is fast ( O(1) )
 - insert_anywhere is fast ( O(1) )
 - mallocs are infrequent
 - memory overhead can be large if frequent deletions and push_backs are performed
 - random access not supported
 - iterators are invalidated by insertions but not deletions.
   iteration however, can be slow if frequent deletions and push_backs are performed.
 - pointers to value_type are never invalidated
 
The implementation follows that of a deque but where deletions are performed
lazily (by marking) and a free list is is used to manage "insert-anywhere"

\note This is at the moment very incomplete.
      I will complete it if there is demand for it.
*/
template <typename T>
class lazy_deque {
 public:
  typedef std::pair<T, size_t> value_type;
 private:
  std::deque<value_type> data;
  std::deque<value_type*> freelist;
  
  // copy not implemented
  lazy_deque<T>& operator=(const lazy_deque<T> &other) { }
  
 public:
  lazy_deque() { }
  
  value_type* push_back(T& dat) {
    data.push_back(lazy_deque_container(dat));
    return &(data.back());
  }
  
  value_type* push_anywhere(T& dat) {
    if (!freelist.empty()) {
     value_type* ret = freelist.front();
     ret->first = dat;
     ret->second = false;
     freelist.pop_front();
     return ret;
    }
    else {
      data.push_back(std::make_pair(dat, false));
      return &(data.back());
    }
  }
  
  void erase(value_type* e) {
    e->second = true;
    freelist.push_back(e);
  }
};

}

#endif
