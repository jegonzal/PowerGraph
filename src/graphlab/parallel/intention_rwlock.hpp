#ifndef INTENTION_LOCK_HPP
#define INTENTION_LOCK_HPP
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#define READERS 0
#define INTENTION 1
#define WRITE 3
#define NONE 4

namespace graphlab {

/**
 Multilevel/Intention Lock implementation
 The lock permits 4 lock types.
  - Read Lock
  - Write Lock
  - Intention Read Lock
  - Intention Write Lock
  
  With the following coexistance rules
 --------------------------------
      |  R  |  W  |  IR  |  IW  |
 -----|-------------------------|
   R  |  Y  |  N  |  Y   |  N   |
 -----|-------------------------|
   W  |  N  |  N  |  N   |  N   |
 -----|-------------------------|
  IR  |  Y  |  N  |  Y   |  Y   |
 -----|-------------------------|
  IW  |  N  |  N  |  Y   |  Y   |
 -----|-------------------------|
 
 This can be used as the upper levels of a multi-level
 lock implementation.
 Essentially, given a lock hierarchy
 
           Parent
            / |  \
           /  |   \
          C1  C2  C3
          
where the parent is an intention lock.
All children C1...C3 can be acquired simultaneously by acquiring a lock on the 
parent (either read or write).

However, to acquire a read lock on just a single child, the caller should 
acquire an intention-read on the parent before acquiring a read on the desired 
child. Similarly for writes.
*/
class intention_lock {
 public:
    intention_lock() {};
    inline void readlock() {
      B.readlock();
    }
    
    
    inline void writelock() {
      B.writelock();
      C.writelock();
    }
    
    
    inline void intention_readlock() {
      C.readlock();
    }
    
    inline void intention_writelock() {
      char nw = num_intention_writes.inc();
      if (nw == 1) {
        B.writelock();
        wr0_acquired.inc();
      }
      else {
        while(wr0_acquired.value == 0);
      }
    }

    inline void rdunlock() {
      B.rdunlock();
    }
    
    
    inline void wrunlock() {
      C.wrunlock();
      B.wrunlock();
    }
    
    inline void intention_rdunlock() {
      C.rdunlock();
    }
    
    inline void intention_wrunlock() {
      char nw = num_intention_writes.inc();
      if (nw == 0) {
        wr0_acquired.dec();
        B.wrunlock();
      }
    }

 private:
    rwlock B;
    rwlock C;
    atomic<uint16_t> num_intention_writes;
    atomic<char> wr0_acquired; /** whether the first thread acquiring an intention write
                                   has acquire lock B. */
};

}

#undef READ_STATE
#undef WRITE_STATE
#undef INTENTION_READ_STATE
#undef INTENTION_WRITE_STATE

#endif
