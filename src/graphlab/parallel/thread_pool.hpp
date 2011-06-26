#ifndef GRAPHLAB_THREAD_POOL_HPP
#define GRAPHLAB_THREAD_POOL_HPP

#include <boost/bind.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/blocking_queue.hpp>

namespace graphlab {

  /**
   * \ingroup util
   * Manages a pool of threads.
   * 
   * The interface is nearly identical to the \ref thread_group. 
   * The key difference is internal behavior. The thread pool preallocates a
   * collection of threads which it keeps asleep. When tasks are issued 
   * through the "launch" function, threads are woken up to perform the
   * tasks. 
   *
   * The thread_pool object performs limited exception forwarding.
   * exception throws within a thread of type const char* will be caught
   * and forwarded to the join() function.
   * If the call to join() is wrapped by a try-catch block, the exception
   * will be caught safely and thread cleanup will be completed properly.
   *
   * If multiple threads are running in the thread-group, the master should
   * test if running_threads() is > 0, and retry the join().
   *
   */
  class thread_pool {
    private:

      thread_group threads;
      blocking_queue<boost::function<void (void)> > spawn_queue;

      // protects the exception queue, and the task counters
      mutex mut;
      conditional event_condition;  // to wake up the joining thread
      std::queue<const char*> exception_queue;
      size_t tasks_inserted;
      size_t tasks_completed;
      bool waiting_on_join; // true if a thread is waiting in join

      // not implemented
      thread_pool& operator=(const thread_pool &thrgrp);
      thread_pool(const thread_pool&);
      
      void wait_for_task();

    public:
      
       /* Initializes a thread pool with nthreads. 
        * If affinity is set, the nthreads will by default stripe across 
        * the available cores on the system. 
        */
      thread_pool(size_t nthreads, bool affinity = false);
  
      /** 
       * Launch a single thread which calls spawn_function. If affinity
       * is set on construction of the thread_pool, the thread handling the
       * function will be locked on to one particular CPU.
       */
      void launch(const boost::function<void (void)> &spawn_function);
  
  
      /** Waits for all threads to become free. const char* exceptions
      thrown by threads are forwarded to the join() function.
      Once this function returns normally, the queue is empty.
      
      Note that this function may not return if producers continually insert
      tasks through launch. 
      */
      void join();
      
      //! Destructor. Cleans up all threads
      ~thread_pool();
  };
  
}
#endif
