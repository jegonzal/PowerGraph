#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/queued_rwlock.hpp>
#include <graphlab/parallel/deferred_rwlock.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/random.hpp>
#include <vector>
#include <queue>
using namespace graphlab;

#define NUM_LOCKS 1000
#define NUM_RAND 500
#define NUM_ITER 10000
deferred_rwlock locks[NUM_LOCKS];
queued_rw_lock queuedlocks[NUM_LOCKS];
rwlock regularlocks[NUM_LOCKS];
atomic<int> readers[NUM_LOCKS];
#define nthreads 32
#define READ_PROP 2
barrier bar1(nthreads);
barrier bar2(nthreads);
void eval_wr(size_t lockid) ;
void eval_rd(size_t lockid);
void release_wr(size_t lockid, deferred_rwlock::request* l) ;
void release_rd(size_t lockid, deferred_rwlock::request* l);
void iterate_released(size_t lockid, deferred_rwlock::request *l, size_t numok) ;

atomic<size_t> numacquired;
atomic<size_t> numreleased;

void eval_wr(size_t lockid) {
  numacquired.inc();
  ASSERT_EQ(readers[lockid].value, 0);
}

void eval_rd(size_t lockid) {
  readers[lockid].inc();
  numacquired.inc();
  readers[lockid].dec();
}

void release_wr(size_t lockid, deferred_rwlock::request* l) {
  deferred_rwlock::request* reqs;
  size_t numok = locks[lockid].wrunlock(reqs);
  numreleased.inc();
  iterate_released(lockid, reqs, numok);
}

void release_rd(size_t lockid, deferred_rwlock::request* l) {
  deferred_rwlock::request* reqs;
  size_t numok = locks[lockid].rdunlock(reqs);
  numreleased.inc();
  iterate_released(lockid, reqs, numok);
}

void iterate_released(size_t lockid, 
		      deferred_rwlock::request *reqs, size_t numok) {
  for (size_t i = 0;i < numok; ++i) {
    if (reqs->lockclass == QUEUED_RW_LOCK_REQUEST_READ) {
      eval_rd(lockid);
      release_rd(lockid, reqs);
    }
    else if (reqs->lockclass == QUEUED_RW_LOCK_REQUEST_WRITE) {
      eval_wr(lockid);
      release_wr(lockid, reqs);
    }
    reqs = (deferred_rwlock::request*)reqs->next;
  }
}


void f(void) {
  std::vector<size_t> randlocks;
  std::vector<bool> randsign;  
  randlocks.resize(NUM_RAND);
  randsign.resize(NUM_RAND);
  deferred_rwlock::request testreqs[NUM_RAND];
    for (size_t j = 0;j < NUM_RAND; ++j) {
      randlocks[j] = graphlab::random::rand_int(NUM_LOCKS - 1);
      randsign[j] = graphlab::random::rand_int(READ_PROP);
    }
for (size_t i = 0;i < NUM_ITER; ++i) {
    bar1.wait();
    std::sort(randlocks.begin(), randlocks.end());
    for (size_t j = 0;j < NUM_RAND; ++j) {
      if (randsign[j] == 0) {
        if (locks[randlocks[j]].writelock(&testreqs[j])) {
          eval_wr(randlocks[j]);
          release_wr(randlocks[j], &testreqs[j]);
        }
      }
      else {
        deferred_rwlock::request *allreleased;
        size_t numok = 
	  locks[randlocks[j]].readlock(&testreqs[j], allreleased);
        iterate_released(randlocks[j], allreleased, numok);
      }
    }
    bar2.wait();
    ASSERT_EQ(numacquired.value, (i+1) * NUM_RAND * nthreads);
    // check that all locks are released
   /* for (size_t j = 0;j < NUM_RAND; ++j) {
      ASSERT_EQ(testreqs[j].s.state.blocked, UNBLOCKED_AND_RETURNED);
    }
    for (size_t j = 0;j < NUM_LOCKS; ++j) {
      ASSERT_EQ(locks[j].get_reader_count(), 0);
      ASSERT_FALSE(locks[j].has_waiters());
      #ifdef DEFERRED_LOCK_DEBUG
      locks[j].nreadpath1.value = 0;
      locks[j].nreadpath2.value = 0;
      locks[j].nwrunlock.value = 0;
      locks[j].nwrpath.value = 0;
      locks[j].nrlock.value = 0;
      locks[j].nwrlock.value = 0;
      locks[j].nreadpath3.value = 0;
      locks[j].nrdunlock.value = 0;
      locks[j].nreadpath4.value = 0;
      locks[j].nreadpath5.value = 0;
      locks[j].nreadpath6.value = 0;
      #endif
    }*/
  }  
}


void f2(void) {
  std::vector<size_t> randlocks;
  std::vector<bool> randsign;  
  randlocks.resize(NUM_RAND);
  randsign.resize(NUM_RAND);
  for (size_t j = 0;j < NUM_RAND; ++j) {
    randlocks[j] = graphlab::random::rand_int(NUM_LOCKS - 1);
    randsign[j] = graphlab::random::rand_int(READ_PROP);
  }
    
  for (size_t i = 0;i < NUM_ITER; ++i) {
    bar1.wait();

    std::sort(randlocks.begin(), randlocks.end());
    for (size_t j = 0;j < NUM_RAND; ++j) {
      if (randsign[j] == 0) {
        regularlocks[randlocks[j]].writelock();
        eval_wr(randlocks[j]);
        regularlocks[randlocks[j]].unlock();
      }
      else {
        regularlocks[randlocks[j]].readlock();
        eval_rd(randlocks[j]);
        regularlocks[randlocks[j]].unlock();
      }
    }
    bar2.wait();
  }  
}




void f3(void) {
  std::vector<size_t> randlocks;
  std::vector<bool> randsign;  
  randlocks.resize(NUM_RAND);
  randsign.resize(NUM_RAND);
  queued_rw_lock::request req[NUM_LOCKS];
  for (size_t j = 0;j < NUM_RAND; ++j) {
    randlocks[j] = graphlab::random::rand_int(NUM_LOCKS - 1);
    randsign[j] = graphlab::random::rand_int(READ_PROP);
  }
  
  for (size_t i = 0;i < NUM_ITER; ++i) {
    bar1.wait();
    std::sort(randlocks.begin(), randlocks.end());
    for (size_t j = 0;j < NUM_RAND; ++j) {
      if (randsign[j] == 0) {
        queuedlocks[randlocks[j]].writelock(&req[randlocks[j]]);
        eval_wr(randlocks[j]);
        queuedlocks[randlocks[j]].wrunlock(&req[randlocks[j]]);
      }
      else {
        queuedlocks[randlocks[j]].readlock(&req[randlocks[j]]);
        eval_rd(randlocks[j]);
        queuedlocks[randlocks[j]].rdunlock(&req[randlocks[j]]);
      }
    }
    bar2.wait();
   /* for (size_t j = 0;j < NUM_LOCKS; ++j) {
      ASSERT_EQ(locks[j].get_reader_count(), 0);
      ASSERT_FALSE(locks[j].has_waiters());
    }*/
  }  
}



int main(int argc, char** argv) {
  deferred_rwlock::request reqs[4];
  
  // basic functionality
  deferred_rwlock lock;
  deferred_rwlock::request* released;
  ASSERT_TRUE(lock.writelock(&reqs[0]));
  ASSERT_EQ(lock.wrunlock(released), 0);
  
  // write reads work
  ASSERT_TRUE(lock.readlock(&reqs[0], released) == 1);
  ASSERT_TRUE(released == &reqs[0]);
  ASSERT_TRUE(lock.readlock(&reqs[1], released) == 1);
  ASSERT_TRUE(released == &reqs[1]);
  ASSERT_TRUE(lock.readlock(&reqs[2], released) == 1);
  ASSERT_TRUE(released == &reqs[2]);
  ASSERT_EQ(lock.rdunlock(released), 0);
  ASSERT_EQ(lock.rdunlock(released), 0);
  ASSERT_EQ(lock.rdunlock(released), 0);
  
  // writes block reads
  ASSERT_TRUE(lock.writelock(&reqs[0]));
  ASSERT_EQ(lock.readlock(&reqs[1], released), 0);
  ASSERT_TRUE(lock.readlock(&reqs[2], released) == 0);
  ASSERT_TRUE(lock.readlock(&reqs[3], released) == 0);
  // unlocking the write will release the reads
  ASSERT_EQ(lock.wrunlock(released), 3);
  ASSERT_TRUE(released != NULL);
  // check to see if I have all of them
  ASSERT_TRUE(released == &reqs[1]);
  released = (deferred_rwlock::request*)released->next;
  ASSERT_TRUE(released == &reqs[2]);
  released = (deferred_rwlock::request*)released->next;
  ASSERT_TRUE(released == &reqs[3]);
  ASSERT_EQ(lock.rdunlock(released), 0);
  ASSERT_EQ(lock.rdunlock(released), 0);
  ASSERT_EQ(lock.rdunlock(released), 0);
  
  
  // reads block writes and writes block writes
  ASSERT_TRUE(lock.readlock(&reqs[0],released) == 1);
  ASSERT_TRUE(released == &reqs[0]);
  ASSERT_FALSE(lock.writelock(&reqs[1]));
  ASSERT_FALSE(lock.writelock(&reqs[2]));
  ASSERT_FALSE(lock.writelock(&reqs[3]));
  // releasing the read only releases one write
  ASSERT_EQ(lock.rdunlock(released), 1);
  ASSERT_TRUE(released == &reqs[1]);
  ASSERT_EQ(lock.wrunlock(released), 1);
  ASSERT_TRUE(released == &reqs[2]);
  ASSERT_EQ(lock.wrunlock(released), 1);
  ASSERT_TRUE(released == &reqs[3]);
  ASSERT_EQ(lock.wrunlock(released), 0);
  ASSERT_TRUE(released == NULL);
  
  ASSERT_FALSE(lock.has_waiters());
  ASSERT_EQ(lock.get_reader_count(), 0);
  // threaded test
    //return 0;
  thread_group group;
  timer ti;
  ti.start();
  for (size_t i = 0;i < nthreads ; ++i) {
    launch_in_new_thread(group, f);
  }
  group.join();
  ASSERT_EQ(numacquired.value, (size_t)(nthreads) * NUM_RAND * NUM_ITER);
  std::cout << (size_t)nthreads * NUM_RAND * NUM_ITER 
	    << " deferred locks acquired and released in " 
	    << ti.current_time() << std::endl;

  thread_group group2;
  ti.start();
  for (size_t i = 0;i < nthreads ; ++i) {
    launch_in_new_thread(group2, f2);
  }
  group2.join();
  std::cout << nthreads * NUM_RAND * NUM_ITER 
	    << " regular locks acquired and released in " 
	    << ti.current_time() << std::endl;
  
  thread_group group3;
  ti.start();
  for (size_t i = 0;i < nthreads ; ++i) {
    launch_in_new_thread(group3, f3);
  }
  group3.join();
  std::cout << (size_t)nthreads * NUM_RAND * NUM_ITER 
	    << " queued locks acquired and released in " 
	    << ti.current_time() << std::endl;
}
