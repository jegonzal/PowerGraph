#ifndef REPLY_INCREMENT_COUNTER_HPP
#define REPLY_INCREMENT_COUNTER_HPP
#include <string>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

namespace graphlab {

class distributed_control;

namespace dc_impl {
/**
A wrapper around a char array. This structure 
is incapable of freeing itself and must be managed externally
*/
struct blob {
  blob(char* c, size_t len):c(c),len(len) { };
  blob():c(NULL), len(0){ };
  
  char *c;
  size_t len;
  
  void save(oarchive& oarc) const {
    oarc << len;
    serialize(oarc, c, len);
  }
 void load(iarchive& iarc) {
    if (c) ::free(c);
    iarc >> len;
    c = (char*) malloc(len);
    deserialize(iarc, c, len);
  }
  
  void free() {
    if (c) {
      ::free(c);
      c = NULL;
      len = 0;
    }
  }
};

/**
Defines a really useful function that performs an atomic
increment of a flag when called. This is useful for waiting
for a reply to a request
*/
struct reply_ret_type{
  atomic<size_t> flag;
  blob val;
  bool usesem;
  semaphore sem;
  reply_ret_type(bool usesem):usesem(usesem) { }
  inline void wait() {
    if (usesem) {
      sem.wait();
    }
    else {
      while(flag.value == 0) sched_yield();
    }
  }
};

}


void reply_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, dc_impl::blob ret);

}

#endif
