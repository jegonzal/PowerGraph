#ifndef GRAPHLAB_PTHREAD_TOOLS_HPP
#define GRAPHLAB_PTHREAD_TOOLS_HPP


#include <cstdlib>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>
#include <signal.h>
#include <sys/time.h>
#include <vector>
#include <list>
#include <queue>
#include <iostream>
#include <boost/function.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/atomic.hpp>

 #undef _POSIX_SPIN_LOCKS
#define _POSIX_SPIN_LOCKS -1


#define __likely__(x)       __builtin_expect((x),1)
#define __unlikely__(x)     __builtin_expect((x),0)




/**
 * \file pthread_tools.hpp A collection of utilities for threading
 */
namespace graphlab {

  /**
   * \class mutex 
   * 
   * Wrapper around pthread's mutex On single core systems mutex
   * should be used.  On multicore systems, spinlock should be used.
   */
  class mutex {
  private:
    // mutable not actually needed
    mutable pthread_mutex_t m_mut;
  public:
    mutex() {
      int error = pthread_mutex_init(&m_mut, NULL);
      ASSERT_TRUE(!error);
    }
    inline void lock() const {
      int error = pthread_mutex_lock( &m_mut  );
      if (error) std::cout << "mutex.lock() error: " << error << std::endl;
      ASSERT_TRUE(!error);
    }
    inline void unlock() const {
      int error = pthread_mutex_unlock( &m_mut );
      ASSERT_TRUE(!error);
    }
    inline bool try_lock() const {
      return pthread_mutex_trylock( &m_mut ) == 0;
    }
    ~mutex(){
      int error = pthread_mutex_destroy( &m_mut );
      ASSERT_TRUE(!error);
    }
    friend class conditional;
  }; // End of Mutex

#if _POSIX_SPIN_LOCKS >= 0
  // We should change this to use a test for posix_spin_locks eventually
  
  // #ifdef __linux__
  /**
   * \class spinlock
   * 
   * Wrapper around pthread's spinlock On single core systems mutex
   * should be used.  On multicore systems, spinlock should be used.
   * If pthread_spinlock is not available, the spinlock will be
   * typedefed to a mutex
   */
  class spinlock {
  private:
    // mutable not actually needed
    mutable pthread_spinlock_t m_spin;
  public:
    spinlock () {
      int error = pthread_spin_init(&m_spin, PTHREAD_PROCESS_PRIVATE);
      ASSERT_TRUE(!error);
    }
  
    inline void lock() const { 
      int error = pthread_spin_lock( &m_spin  );
      ASSERT_TRUE(!error);
    }
    inline void unlock() const {
      int error = pthread_spin_unlock( &m_spin );
      ASSERT_TRUE(!error);
    }
    inline bool try_lock() const {
      return pthread_spin_trylock( &m_spin ) == 0;
    }
    ~spinlock(){
      int error = pthread_spin_destroy( &m_spin );
      ASSERT_TRUE(!error);
    }
    friend class conditional;
  }; // End of spinlock
#define SPINLOCK_SUPPORTED 1
#else
  //! if spinlock not supported, it is typedef it to a mutex.
  typedef mutex spinlock;
#define SPINLOCK_SUPPORTED 0
#endif

  
  
  class simple_spinlock {
  private:
    // mutable not actually needed
    mutable volatile char spinner;
  public:
    simple_spinlock () {
      spinner = 0;
    }
  
    inline void lock() const { 
      while(spinner == 1 || __sync_lock_test_and_set(&spinner, 1));
    }
    inline void unlock() const {
      __sync_synchronize();
      spinner = 0;
    }
    inline bool try_lock() const {
      return (__sync_lock_test_and_set(&spinner, 1) == 0);
    }
    ~simple_spinlock(){
      ASSERT_TRUE(spinner == 0);
    }
  };
  

  /**
   * \class conditional
   * Wrapper around pthread's condition variable
   */
  class conditional {
  private:
    mutable pthread_cond_t  m_cond;
  public:
    conditional() {
      int error = pthread_cond_init(&m_cond, NULL);
      ASSERT_TRUE(!error);
    }
    inline void wait(const mutex& mut) const {
      int error = pthread_cond_wait(&m_cond, &mut.m_mut);
      ASSERT_TRUE(!error);
    }
    inline int timedwait(const mutex& mut, int sec) const {
      struct timespec timeout;
      struct timeval tv;
      struct timezone tz;
      gettimeofday(&tv, &tz);
      timeout.tv_nsec = 0;
      timeout.tv_sec = tv.tv_sec + sec;
      return pthread_cond_timedwait(&m_cond, &mut.m_mut, &timeout);
    }
    inline int timedwait_ns(const mutex& mut, int ns) const {
      struct timespec timeout;
      struct timeval tv;
      struct timezone tz;
      gettimeofday(&tv, &tz);
      timeout.tv_nsec = (tv.tv_usec * 1000 + ns) % 1000000000;
      timeout.tv_sec = tv.tv_sec + (tv.tv_usec * 1000 + ns >= 1000000000);
      return pthread_cond_timedwait(&m_cond, &mut.m_mut, &timeout);
    }

    inline void signal() const {
      int error = pthread_cond_signal(&m_cond);
      ASSERT_TRUE(!error);
    }
    inline void broadcast() const {
      int error = pthread_cond_broadcast(&m_cond);
      ASSERT_TRUE(!error);
    }
    ~conditional() {
      int error = pthread_cond_destroy(&m_cond);
      ASSERT_TRUE(!error);
    }
  }; // End conditional

  /**
   * \class semaphore
   * Wrapper around pthread's semaphore
   */
#ifdef __APPLE__
  class semaphore {
  private:
    conditional cond;
    mutex mut;
    mutable volatile size_t semvalue;
    mutable volatile size_t waitercount;
  public:
    semaphore() {
      semvalue = 0;
      waitercount = 0;
    }
    inline void post() const {
      mut.lock();
      if (waitercount > 0) {
        cond.signal();
      }
      semvalue++;
      mut.unlock();
    }
    inline void wait() const {
      mut.lock();
      waitercount++;
      while (semvalue == 0) {
        cond.wait(mut);
      }
      waitercount--;
      semvalue--;
      mut.unlock();
    }
    ~semaphore() {
      ASSERT_TRUE(waitercount == 0);
      ASSERT_TRUE(semvalue == 0);
    }
  }; // End semaphore
#else
  class semaphore {
  private:
    mutable sem_t  m_sem;
  public:
    semaphore() {
      int error = sem_init(&m_sem, 0,0);
      ASSERT_TRUE(!error);
    }
    inline void post() const {
      int error = sem_post(&m_sem);
      ASSERT_TRUE(!error);
    }
    inline void wait() const {
      int error = sem_wait(&m_sem);
      ASSERT_TRUE(!error);
    }
    ~semaphore() {
      int error = sem_destroy(&m_sem);
      ASSERT_TRUE(!error);
    }
  }; // End semaphore
#endif
  

#define atomic_xadd(P, V) __sync_fetch_and_add((P), (V))
#define cmpxchg(P, O, N) __sync_val_compare_and_swap((P), (O), (N))
#define atomic_inc(P) __sync_add_and_fetch((P), 1)

  /**
   * \class spinrwlock
   * rwlock built around "spinning"
   * source adapted from http://locklessinc.com/articles/locks/
   * "Scalable Reader-Writer Synchronization for Shared-Memory Multiprocessors"
   * John Mellor-Crummey and Michael Scott
   */
  class spinrwlock {

    union rwticket {
      unsigned u;
      unsigned short us;
      __extension__ struct {
        unsigned char write;
        unsigned char read;
        unsigned char users;
      } s;
    };
    mutable bool writing;
    mutable volatile rwticket l;
  public:
    spinrwlock() {
      memset(const_cast<rwticket*>(&l), 0, sizeof(rwticket));
    }
    inline void writelock() const {
      unsigned me = atomic_xadd(&l.u, (1<<16));
      unsigned char val = me >> 16;
    
      while (val != l.s.write) sched_yield();
      writing = true;
    }

    inline void wrunlock() const{
      rwticket t = *const_cast<rwticket*>(&l);

      t.s.write++;
      t.s.read++;
    
      *(volatile unsigned short *) (&l) = t.us;
      writing = false;
      __asm("mfence");
    }

    inline void readlock() const {
      unsigned me = atomic_xadd(&l.u, (1<<16));
      unsigned char val = me >> 16;
    
      while (val != l.s.read) sched_yield();
      l.s.read++;
    }

    inline void rdunlock() const {
      atomic_inc(&l.s.write);
    }
  
    inline void unlock() const {
      if (!writing) rdunlock();
      else wrunlock();
    }
  };

#undef atomic_xadd
#undef cmpxchg
#undef atomic_inc


  /**
   * \class rwlock
   * Wrapper around pthread's rwlock
   */
  class rwlock {
  private:
    mutable pthread_rwlock_t m_rwlock;
  public:
    rwlock() {
      int error = pthread_rwlock_init(&m_rwlock, NULL);
      ASSERT_TRUE(!error);
    }
    ~rwlock() {
      int error = pthread_rwlock_destroy(&m_rwlock);
      ASSERT_TRUE(!error);
    }
    inline void readlock() const {
      pthread_rwlock_rdlock(&m_rwlock);
      //ASSERT_TRUE(!error);
    }
    inline void writelock() const {
      pthread_rwlock_wrlock(&m_rwlock);
      //ASSERT_TRUE(!error);
    }
    inline void unlock() const {
      pthread_rwlock_unlock(&m_rwlock);
      //ASSERT_TRUE(!error);
    }
    inline void rdunlock() const {
      unlock();
    }
    inline void wrunlock() const {
      unlock();
    }
  }; // End rwlock

  /**
   * \class barrier
   * Wrapper around pthread's barrier
   */
#ifdef __linux__
  /**
   * \class barrier
   * Wrapper around pthread's barrier
   */
  class barrier {
  private:
    mutable pthread_barrier_t m_barrier;
  public:
    barrier(size_t numthreads) { pthread_barrier_init(&m_barrier, NULL, numthreads); }
    ~barrier() { pthread_barrier_destroy(&m_barrier); }
    inline void wait() const { pthread_barrier_wait(&m_barrier); }
  };

#else
  /**
   * \class barrier
   * Wrapper around pthread's barrier
   */
  class barrier {
  private:
    mutex m;
    int needed;
    int called;
    conditional c;
    
    bool barrier_sense;
    bool barrier_release;
    // we need the following to protect against spurious wakeups
  
  public:
    
    barrier(size_t numthreads) {
      needed = numthreads;
      called = 0;
      barrier_sense = false;
      barrier_release = true;
    }
    
    ~barrier() {}
    
    
    inline void wait() {
      m.lock();
      // set waiting;
      called++;
      bool listening_on = barrier_sense;
      
      if (called == needed) {
        // if I have reached the required limit, wait up. Set waiting
        // to 0 to make sure everyone wakes up

        called = 0;
        barrier_release = barrier_sense;
        barrier_sense = !barrier_sense;
        // clear all waiting
        c.broadcast();
      }
      else {
        // while no one has broadcasted, sleep
        while(barrier_release != listening_on) c.wait(m);
      }
      m.unlock();
    }
  };
#endif



  inline void prefetch_range(void *addr, size_t len) {
    char *cp;
    char *end = (char*)(addr) + len;

    for (cp = (char*)(addr); cp < end; cp += 64) __builtin_prefetch(cp, 0); 
  }
  inline void prefetch_range_write(void *addr, size_t len) {
    char *cp;
    char *end = (char*)(addr) + len;

    for (cp = (char*)(addr); cp < end; cp += 64) __builtin_prefetch(cp, 1);
  }









  /**
   * \class thread 
   * A collection of routines for starting and managing threads.
   * 
   */
  class thread {
  public:

    /**
     * This class contains the data unique to each thread. All threads
     * are gauranteed to have an associated graphlab thread_specific
     * data. The thread object is copyable. 
     */  
    class tls_data {
    public:
      inline tls_data(size_t thread_id) : thread_id_(thread_id) { }
      inline size_t thread_id() { return thread_id_; }
    private:
      size_t thread_id_;
    }; // end of thread specific data



    /// Static helper routines
    // ===============================================================

    /**
     * Get the thread specific data associated with this thread
     */
    static tls_data& get_tls_data();
      
    /** Get the id of the calling thread.  This will typically be the
        index in the thread group. Between 0 to ncpus. */
    static inline size_t thread_id() { return get_tls_data().thread_id(); }
    

    
    /**
     * This static method joins the invoking thread with the other
     * thread object.  This thread will not return from the join
     * routine until the other thread complets it run.
     */
    static void join(thread& other);
    
    // Called just before thread exits. Can be used
    // to do special cleanup... (need for Java JNI)
    static void thread_destroy_callback();
    static void set_thread_destroy_callback(void (*callback)());

      
    /**
     * Return the number processing units (individual cores) on this
     * system
     */
    static size_t cpu_count();


  private:
    
    struct invoke_args{
      size_t m_thread_id;
      boost::function<void(void)> spawn_routine;   
      invoke_args(size_t m_thread_id, const boost::function<void(void)> &spawn_routine)
          : m_thread_id(m_thread_id), spawn_routine(spawn_routine) { };
    };
    
    //! Little helper function used to launch threads
    static void* invoke(void *_args);   
  
  public:
    
    /**
     * Creates a thread with a user-defined associated thread ID
     */
    inline thread(size_t thread_id = 0) : 
      m_stack_size(0), 
      m_p_thread(0),
      m_thread_id(thread_id),
      thread_started(false){
      // Calculate the stack size in in bytes;
      const int BYTES_PER_MB = 1048576; 
      const int DEFAULT_SIZE_IN_MB = 8;
      m_stack_size = DEFAULT_SIZE_IN_MB * BYTES_PER_MB;
    }

    /**
     * execute this function to spawn a new thread running spawn_function
     * routine 
     */
    void launch(const boost::function<void (void)> &spawn_routine);

    /**
     * Same as launch() except that you can specify a CPU on which to
     * run the thread.  This only currently supported in Linux and if
     * invoked on a non Linux based system this will be equivalent to
     * start().
     */
     void launch(const boost::function<void (void)> &spawn_routine, size_t cpu_id);


    /**
     * Join the calling thread with this thread.
     */
    inline void join() {
      if(this == NULL) {
        std::cout << "Failure on join()" << std::endl;
        exit(EXIT_FAILURE);
      }
      join(*this);
    }

    inline bool active() const {
      return thread_started;
    }
    
    inline ~thread() {  }

    inline pthread_t pthreadid() {
      return m_p_thread;
    }
  private:
    
    
    //! The size of the internal stack for this thread
    size_t m_stack_size;
    
    //! The internal pthread object
    pthread_t m_p_thread;
    
    //! the threads id
    size_t m_thread_id;
    
    bool thread_started;
  }; // End of class thread

  



  /**
   * \class thread_group Manages a collection of threads
   * This class is not copyable
   */
  class thread_group {
   private:
    size_t m_thread_counter;
    size_t threads_running;
    mutex mut;
    conditional cond;
    std::queue<std::pair<pthread_t, const char*> > joinqueue;
    // not implemented
    thread_group& operator=(const thread_group &thrgrp);
    thread_group(const thread_group&);
    static void invoke(boost::function<void (void)> spawn_function, thread_group *group);
   public:
    /** 
     * Initializes a thread group. 
     */
    thread_group() : m_thread_counter(0), threads_running(0) { }

    /** 
     * Launch a single thread which calls spawn_function No CPU affinity is
     * set so which core it runs on is up to the OS Scheduler
     */
    void launch(const boost::function<void (void)> &spawn_function);

    /**
     * Launch a single thread which calls spawn_function Also sets CPU
     *  Affinity
     */
    void launch(const boost::function<void (void)> &spawn_function, size_t cpu_id);

    //! Waits for all threads to complete execution
    void join();
    
    inline size_t running_threads() {
      return threads_running;
    }
    //! Destructor. Waits for all threads to complete execution
    inline ~thread_group(){ join(); }

  }; // End of thread group


  /// Runs f in a new thread. convenience function for creating a new thread quickly.
  inline thread launch_in_new_thread(const boost::function<void (void)> &f, 
                               size_t cpuid = size_t(-1)) {
    thread thr;
    if (cpuid != size_t(-1)) thr.launch(f, cpuid);
    else thr.launch(f);
    return thr;
  }



}; // End Namespace

#endif
