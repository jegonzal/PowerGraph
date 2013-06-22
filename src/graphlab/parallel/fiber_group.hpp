#ifndef GRAPHLAB_FIBER_GROUP_HPP
#define GRAPHLAB_FIBER_GROUP_HPP
#include <graphlab/parallel/fiber_control.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
namespace graphlab {

/**
 * Defines a group of fibers. Analogous to the thread_group, but is meant
 * to run only little user-mode threads. It is important that fibers never 
 * block, since there is no way to context switch out from a blocked fiber.
 * The fiber_group uses the fiber_control singleton instance to manage its
 * fibers.
 */
class fiber_group {
 public:
  typedef fiber_control::affinity_type affinity_type;

 private:
  size_t stacksize;
  affinity_type affinity;
  atomic<size_t> threads_running;
  mutex join_lock;
  // to be triggered once the threads_running counter becomes 0
  conditional join_cond; 
  // set to true if someone is waiting on a join()
  bool join_waiting;

  inline void increment_running_counter() {
    threads_running.inc();
  }

  inline void decrement_running_counter() {
    // now, a bit of care is needed here
    size_t r = threads_running.dec();
    if (r == 0) {
      join_lock.lock();
      if (join_waiting) {
        join_cond.signal();
      }
      join_lock.unlock();
    }
  }

  // wraps the call so that we can do the appropriate termination
  static void invoke(const boost::function<void (void)>& spawn_function, 
                     fiber_group* group);

 public:


  fiber_group(size_t stacksize = 8192, 
              affinity_type affinity = fiber_control::all_affinity()) : 
      stacksize(stacksize), 
      affinity(affinity),
      join_waiting(false) { }


  /**
   * Sets the stacksize of each fiber.
   * Only takes effect for threads launched after this.
   */
  inline void set_stacksize(size_t new_stacksize) {
    stacksize = new_stacksize;
  }


  /**
   * Sets the affinity for each fiber.
   * Only takes effect for threads launched after this.
   */
  inline void set_affinity(affinity_type new_affinity) {
    affinity = new_affinity;
  }

  /**
   * Launch a single thread which calls spawn_function.
   */
  void launch(const boost::function<void (void)> &spawn_function);

  /** Waits for all threads to complete execution. const char* exceptions
   *  thrown by threads are forwarded to the join() function.
   */
  void join();

  /// Returns the number of running threads.
  inline size_t running_threads() {
    return threads_running;
  }
  //
  //! Destructor. Waits for all threads to complete execution
  inline ~fiber_group(){ join(); }

};

} // namespace graphlab 
#endif
