#ifndef LOCK_FREE_POOL_HPP
#define LOCK_FREE_POOL_HPP
#include <stdint.h>
#include <vector>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/atomic.hpp>

namespace graphlab {
template <typename T, typename index_type = uint32_t>
class lock_free_pool{
 private:
  std::vector<T> data;
  
  // freelist[i] points to the next free list element
  // if freelist[i] == index_type(-1), then it is the last element
  // allocated entries are set to index_type(0), though
  // note that there is no way to disambiguate between allocated
  // and non-allocated entries by simply looking at the freelist
  std::vector<index_type> freelist;
  volatile index_type freelisthead;
 public:
  lock_free_pool(size_t poolsize = 0) {
   reset_pool(poolsize);
  }
  
  
  void reset_pool(size_t poolsize) {
    if (poolsize == 0) {
     data.clear();
     freelist.clear();
    }
    else {
      data.resize(poolsize);
      freelist.resize(poolsize);
      for (index_type i = 0;i < freelist.size(); ++i) {
        freelist[i] = i + 1;
      }
      freelist[freelist.size() - 1] = index_type(-1);
    }
    freelisthead = 0;
  }
  
  std::vector<T>& unsafe_get_pool_ref() {
    return data;
  }
  
  T* alloc() {
    // I need to atomically advance freelisthead to the freelist[head]
    index_type oldhead;
    index_type newhead;
    do {
      oldhead = freelisthead;
      if (oldhead == index_type(-1)) return NULL;
      newhead = freelist[oldhead];
    } while(!atomic_compare_and_swap(freelisthead, oldhead, newhead));
    freelist[oldhead] = index_type(0);
    return &(data[oldhead]);
  }
  
  void free(T* p) {
    // sanity check
    index_type cur = index_type(p - &(data[0]));
    ASSERT_EQ(freelist[cur], index_type(0));
    // prepare for free list insertion
    // I need to atomically set freelisthead == cur
    // and freelist[cur] = freelisthead
    index_type oldhead;
    do{
      oldhead = freelisthead;
      freelist[cur] = oldhead;
      // now try to atomically move freelisthead
    } while(!atomic_compare_and_swap(freelisthead, oldhead, cur));
  }
};

}
#endif
