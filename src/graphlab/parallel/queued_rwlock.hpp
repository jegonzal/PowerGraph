#ifndef DEFERRED_QUEUED_RWLOCK_HPP
#define DEFERRED_QUEUED_RWLOCK_HPP


namespace graphlab {
  
  
#define QUEUED_RW_LOCK_REQUEST_READ 0
#define QUEUED_RW_LOCK_REQUEST_WRITE 1
#define QUEUED_RW_LOCK_REQUEST_NONE 2

/**
 * Fair rw-lock with local-only spinning implemented and
 * modified from 
 * Scalable Reader-Writer Synchronization for Shared-Memory Multiprocessors.
 * John M. Mellor-Crummey and Michael L. Scott
 */
class queued_rw_lock{
 public:
  
  union state_union {
    volatile uint32_t stateu;
    struct {
      volatile uint16_t successor_class;
      volatile bool blocked;
    } state;
  };
  
  struct request{
    void* id;  
    volatile request* volatile next;
    volatile state_union s;
    volatile char lockclass;
  };
 private:
  request* volatile tail;
  atomic<size_t> reader_count;
  request* volatile next_writer;
 public:
  queued_rw_lock(): tail(NULL), reader_count(0), next_writer(NULL) { }
  
  inline void writelock(request *I) {
    I->lockclass = QUEUED_RW_LOCK_REQUEST_WRITE;
    I->next = NULL;
    I->s.stateu = 0;
    I->s.state.blocked = true;
    I->s.state.successor_class = QUEUED_RW_LOCK_REQUEST_NONE;
    __sync_synchronize();
    request* predecessor = __sync_lock_test_and_set(&tail, I);

    if (predecessor == NULL) {
      next_writer = I;
      __sync_synchronize();
      if (reader_count.value == 0 && __sync_lock_test_and_set(&next_writer, (request*)NULL) == I) {
        I->s.state.blocked = false;
      }
    }
    else {
      predecessor->s.state.successor_class = QUEUED_RW_LOCK_REQUEST_WRITE;
    __sync_synchronize();      
      predecessor->next = I;
    }
    // while I->blocked. continue
    volatile state_union& is = I->s;
    while (is.state.blocked) sched_yield();
  }

  inline void wrunlock(request *I) {
    __sync_synchronize(); 

    if (I->next != NULL || !__sync_bool_compare_and_swap(&tail, I, (request*)NULL)) {
      // wait
      while(I->next == NULL) sched_yield();
     __sync_synchronize(); 
   
      if (I->next->lockclass == QUEUED_RW_LOCK_REQUEST_READ) {
        reader_count.inc();
      }
      I->next->s.state.blocked = false;
    }
  }

  inline void readlock(request *I)  {
    I->lockclass =QUEUED_RW_LOCK_REQUEST_READ;
    I->next = NULL;
    I->s.stateu = 0;
    I->s.state.successor_class = QUEUED_RW_LOCK_REQUEST_NONE;
    I->s.state.blocked = true;
    __sync_synchronize(); 
    request* predecessor = __sync_lock_test_and_set(&tail, I);

    if (predecessor == NULL) {
      reader_count.inc();
      I->s.state.blocked = false;
    }
    else {
      
      state_union tempold, tempnew;
      tempold.state.blocked = true;
      tempold.state.successor_class = QUEUED_RW_LOCK_REQUEST_NONE;
      tempnew.state.blocked = true;
      tempnew.state.successor_class = QUEUED_RW_LOCK_REQUEST_READ;
      __sync_synchronize();
      if (predecessor->lockclass == QUEUED_RW_LOCK_REQUEST_WRITE || atomic_compare_and_swap(predecessor->s.stateu, tempold.stateu, tempnew.stateu)) {
        
        predecessor->next = I;
        // wait
        __sync_synchronize(); 
        volatile state_union& is = I->s;
        while(is.state.blocked) sched_yield();
      }
      else {
        reader_count.inc();
        predecessor->next = I;
        __sync_synchronize();
        I->s.state.blocked = false;
      }
    }
    __sync_synchronize();
    if (I->s.state.successor_class == QUEUED_RW_LOCK_REQUEST_READ) {
      
      // wait
      while(I->next == NULL) sched_yield();
      reader_count.inc();
      I->next->s.state.blocked = false;
    }
  }

  inline void rdunlock(request *I)  {
    __sync_synchronize();
    if (I->next != NULL || !__sync_bool_compare_and_swap(&tail, I, (request*)NULL)) {
      while(I->next == NULL) sched_yield();
      if (I->s.state.successor_class == QUEUED_RW_LOCK_REQUEST_WRITE) {
        next_writer = (request*)(I->next);
        __sync_synchronize();
      }
    }
    if (reader_count.dec() == 0) {
      __sync_synchronize(); 
      request * w = __sync_lock_test_and_set(&next_writer, (request*)NULL);
      if (w != NULL) {
        w->s.state.blocked = false;
      }
    }
  }
};
















class deferred_rw_lock{
 public:
  
  union state_union {
    uint32_t stateu;
    struct {
      uint16_t successor_class;
      uint8_t blocked;  
      uint8_t returned;
    } state;
    
  struct {
      uint16_t successor_class;
      uint16_t blocked_and_returned;
    } state2;
  };
  
  struct request{
    void* id;  
    volatile request* volatile next;
    volatile state_union s;
    char lockclass;
  };
 private:
  request* tail;
  atomic<size_t> reader_count;
  request* volatile next_writer;

 public:
  //atomic<size_t> nreadpath1, nreadpath2;
  deferred_rw_lock(): tail(NULL), reader_count(0), next_writer(NULL) { }

  // debugging purposes only
  size_t get_reader_count() {
    return reader_count.value;
  }
  
  // debugging purposes only
  bool has_waiters() {
    return tail != NULL;
  }
  
  inline bool writelock(request *I) {
    I->lockclass = QUEUED_RW_LOCK_REQUEST_WRITE;
    I->next = NULL;
    I->s.stateu = 0;
    I->s.state.blocked = true;
    I->s.state.returned = false;
    I->s.state.successor_class = QUEUED_RW_LOCK_REQUEST_NONE;

    request* predecessor = fetch_and_store(tail, I);

    if (predecessor == NULL) {
      next_writer = I;
      __sync_synchronize();
      if (reader_count.value == 0 &&  __sync_lock_test_and_set(&next_writer, (request*)NULL) == I) {
        I->s.state.blocked = false;
      }
    }
    else {
      predecessor->s.state.successor_class = QUEUED_RW_LOCK_REQUEST_WRITE;
    __sync_synchronize();      
      predecessor->next = I;
    }
    // if not blocked, set returned and return
    volatile state_union& is = I->s;
    state_union unblocked_not_returned, unblocked_returned;
    unblocked_not_returned.state.blocked = false;
    unblocked_not_returned.state.returned = false;
    unblocked_returned.state.blocked = false;
    unblocked_returned.state.returned = true;
    return (atomic_compare_and_swap(is.state2.blocked_and_returned, 
                                unblocked_not_returned.state2.blocked_and_returned, 
                                unblocked_returned.state2.blocked_and_returned));
  }


  inline bool resume_wrlock(request* I) {
    // to be called when blocked == false
    volatile state_union& is = I->s;
    state_union unblocked_not_returned, unblocked_returned;
    unblocked_not_returned.state.blocked = false;
    unblocked_not_returned.state.returned = false;
    unblocked_returned.state.blocked = false;
    unblocked_returned.state.returned = true;
    if (atomic_compare_and_swap(is.state2.blocked_and_returned, 
                                unblocked_not_returned.state2.blocked_and_returned, 
                                unblocked_returned.state2.blocked_and_returned)) return true;
    else return false;  // the lock has already been returned
  }
  
  inline size_t wrunlock(request *I, request* &released) {
    //nwrite.inc();
    released = NULL;
    size_t numreleased  = 0;
    if (I->next != NULL || !atomic_compare_and_swap(tail, I, (request*)NULL)) {
      // wait
      while(I->next == NULL) sched_yield();
      // unlock the next sequence
      
      
      if (I->next->lockclass == QUEUED_RW_LOCK_REQUEST_READ) {
        I->next->s.state.blocked = false;
        numreleased = resume_rdlock((request*)(I->next), released);
      }
      else if (I->next->lockclass == QUEUED_RW_LOCK_REQUEST_WRITE) {
        I->next->s.state.blocked = false;
        if (resume_wrlock((request*)(I->next))) {
          released = (request*)(I->next);
          numreleased  = 1;
        }
      }
    }
    return numreleased ;
  }

  inline size_t readlock(request *I, request* &allreleases)  {
    allreleases = NULL;
    size_t numreleased = 0;
    I->lockclass =QUEUED_RW_LOCK_REQUEST_READ;
    I->next = NULL;
    I->s.stateu = 0;
    I->s.state.successor_class = QUEUED_RW_LOCK_REQUEST_NONE;
    I->s.state.returned = false;
    I->s.state.blocked = true;
    
    request* predecessor = fetch_and_store(tail, I);
    if (predecessor == NULL) {
      reader_count.inc();
      I->s.state.returned = true;
      __sync_synchronize();
      I->s.state.blocked = false;
      allreleases = I;
      numreleased = 1;
    }
    else {
      
      state_union tempold, tempnew;
      tempold.state.blocked = true;
      tempold.state.returned = false;
      tempold.state.successor_class = QUEUED_RW_LOCK_REQUEST_NONE;
      tempnew.state.blocked = true;
      tempnew.state.returned = false;
      tempnew.state.successor_class = QUEUED_RW_LOCK_REQUEST_READ;
      if (predecessor->lockclass == QUEUED_RW_LOCK_REQUEST_WRITE ||  atomic_compare_and_swap(predecessor->s.stateu, tempold.stateu, tempnew.stateu)) {
        
        predecessor->next = I;
        // wait
        __sync_synchronize(); 
        volatile state_union& is = I->s;
        state_union unblocked_not_returned, unblocked_returned;
        unblocked_not_returned.state.blocked = false;
        unblocked_not_returned.state.returned = false;
        unblocked_returned.state.blocked = false;
        unblocked_returned.state.returned = true;
        if (atomic_compare_and_swap(is.state2.blocked_and_returned, 
                                    unblocked_not_returned.state2.blocked_and_returned, 
                                    unblocked_returned.state2.blocked_and_returned)) {
          allreleases = I;
          numreleased = 1;
          reader_count.inc(); 
        }

      }
      else {
        reader_count.inc();
        predecessor->next = I;
        I->s.state.returned = true;
        __sync_synchronize();
        I->s.state.blocked = false;
        allreleases = I;
        numreleased = 1;
      }
    }
    
    
    if (numreleased > 0 && I->s.state.successor_class == QUEUED_RW_LOCK_REQUEST_READ) {
      // wait
      while(I->next == NULL) sched_yield();
      I->next->s.state.blocked = false;
      request* unused;
      numreleased += resume_rdlock((request*)(I->next), unused);
    }
    

    return numreleased;
  }


  inline size_t resume_rdlock(request *I, request* &released) {
    released = NULL;
    size_t numreleased = 0;
    volatile state_union& is = I->s;
    state_union unblocked_not_returned, unblocked_returned;
    unblocked_not_returned.state.blocked = false;
    unblocked_not_returned.state.returned = false;
    unblocked_returned.state.blocked = false;
    unblocked_returned.state.returned = true;
    if (atomic_compare_and_swap(is.state2.blocked_and_returned, 
                                unblocked_not_returned.state2.blocked_and_returned, 
                                unblocked_returned.state2.blocked_and_returned)) {
      reader_count.inc();
      released = I;
      ++numreleased ;
    }
    
    if (I->s.state.successor_class == QUEUED_RW_LOCK_REQUEST_READ) {
      while(I->next == NULL) sched_yield();
      I->next->s.state.blocked = false;
      request* unused;
      numreleased += resume_rdlock((request*)(I->next), unused);
      if (released == NULL) released = unused;
    }
    return numreleased;
  }

  inline size_t rdunlock(request *I, request* &released)  {
    //nread.inc();
    released = NULL;
    size_t numreleased = 0;
    if (I->next != NULL || !atomic_compare_and_swap(tail, I, (request*)NULL)) {
      while(I->next == NULL) sched_yield();
//      nreadpath1.inc();
      if (I->s.state.successor_class == QUEUED_RW_LOCK_REQUEST_WRITE) {
        next_writer = (request*)(I->next);
        __sync_synchronize();
      }
    }
    if (reader_count.dec() == 0) {
      //nreadpath2.inc();
      request * w = __sync_lock_test_and_set(&next_writer, (request*)NULL);
      if (w != NULL) {
        w->s.state.blocked = false;
        if (resume_wrlock(w)) {
          released = w;
          numreleased = 1;
        }
      }
    }
    return numreleased;
  }
};
}
#endif
