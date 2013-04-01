#ifndef GRAPHLAB_FIBER_HPP
#define GRAPHLAB_FIBER_HPP
#include <cstdlib>
#include <boost/context/all.hpp>
#include <boost/function.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
namespace graphlab {

/// a very simple user-mode threading
class fiber_group {
 public:
  struct fiber {
    simple_spinlock lock;
    fiber_group* parent;
    boost::context::fcontext_t* context;
    void* stack;
    size_t id;
    void* fls; // fiber local storage
    fiber* next;
    intptr_t initial_trampoline_args;
    pthread_mutex_t* deschedule_lock; // if descheduled is set, we will
                                      // atomically deschedule and unlock
                                      // this mutex
    bool descheduled; // flag. set if this fiber is to be descheduled.
                      // This is a temporary flag, and is only used to notify
                      // the context switch to deschedule this thread.
                      // lock must be acquired for this to be modified

    bool terminate;   // flag. set if this fiber is to be destroyed.
                      // This is a temporary flag, and is only used to notify
                      // the context switch to destroy this thread.

    bool scheduleable;     // Managed by the queue management routines.
                      // Set if the fiber is inside the scheduling queue
                      // or is running in a thread.
                      // lock must be acquired for this to be modified.
  };


 private:
  size_t nworkers;
  size_t stacksize;
  atomic<size_t> fiber_id_counter;
  atomic<size_t> fibers_active;
  mutex join_lock;
  conditional join_cond;

  bool stop_workers;

 // The scheduler is a single queue
  mutex active_lock;
  conditional active_cond;
  // used as a sentinel for the actual queues
  // first element is always head,
  // tail points to the last element
  fiber active_head;
  fiber* active_tail;
  size_t nactive;

  thread_group workers;


  // locks must be acquired outside the call
  void active_queue_insert(fiber* value);
  fiber* active_queue_remove();

  // a thread local storage for the worker to point to a fiber
  static bool tls_created;
  struct tls {
    fiber_group* parent;
    fiber* prev_fiber; // the fiber we context switch from
    fiber* cur_fiber; // the fiber we are context switching to
    fiber* garbage;
    size_t worker_id;
    boost::context::fcontext_t base_context;
  };
  static pthread_key_t tlskey; // points to the tls structure above
  static void tls_deleter(void* tls);

  static void create_tls_ptr();
  static tls* get_tls_ptr();
  static fiber* get_active_fiber();
  static void set_active_fiber(fiber*);

  void worker_init(size_t workerid);

  void reschedule_fiber(fiber* pfib);
  void yield_to(fiber* next_fib);
  static void trampoline(intptr_t _args);

  void (*flsdeleter)(void*);
 public:
  /// creates a group of fibers using a certain number of worker threads.
  fiber_group(size_t nworkers, size_t stacksize);

  ~fiber_group();

  /** the basic launch function
   * Returns a fiber ID. IDs are not sequential.
   * \note The ID is really a pointer to a fiber_group::fiber object.
   */
  size_t launch(boost::function<void (void)> fn);

  /**
   * Waits for all functions to join
   */
  void join();


  /**
   * Returns the number of threads that have yet to join
   */
  inline size_t num_threads() {
    return fibers_active.value;
  }

  /**
   * Returns the total number threads ever created
   */
  inline size_t total_threads_created() {
    return fiber_id_counter.value;
  }
  /**
   * Sets the TLS deletion function. The deletion function will be called
   * on every non-NULL TLS value.
   */
  void set_tls_deleter(void (*deleter)(void*));
  /**
   * Gets the TLS value. Defaults to NULL.
   * Note that this function will only work within a fiber.
   */
  static void* get_tls();
  /**
   * Sets the TLS value.
   * Note that this function will only work within a fiber.
   * If the value is not NULL, and the deletion function is set by
   * set_tls_deleter(), the deleter will be called on the value on
   * fiber termination.
   */
  static void set_tls(void* value);
  /**
   * Kills the current fiber.
   * Note that this function will only work within a fiber.
   */
  static void exit();
  /**
   * Yields to another fiber.
   * Note that this function will only work within a fiber.
   */
  static void yield();
  /**
   * Returns the current fiber handle.
   * Note that fiber handles are not sequential, and are really a
   * pointer to an internal datastructure.
   * The fiber handle will never be 0.
   * This function will only work within a fiber.
   */
  static size_t get_tid();

  /**
   * Returns the worker managing the current fiber.
   * This function will only work within a fiber.
   * Worker IDs are sequential.
   */
  static size_t get_worker_id();


  /**
   * Atomically deschedules the current thread and unlocks the mutex.
   *
   * deschedule_self() and schedule_tid() must be managed carefully
   * to avoid race conditions. i.e. schedule_tid() happending before
   * deschedule_self().
   *
   * To support this correctly, the descheduling must be paired together
   * with a mutex.
   *
   * For instance, to use this to implement a promise.
   * \code
   * // descheduling fiber
   * pthread_mutex_lock(&lock);
   * if ( ... promise not ready ...) {
   *   deschedule_self(&lock);
   * } else {
   *   pthread_mutex_unlock(&lock);
   * }
   *   ... use the promise ...
   * \endcode
   *
   *
   * The promise execution thread then must do the following
   * \code
   * ... tid contains the fiber ID to wake when promise is done
   * ... compute promise...
   * pthread_mutex_lock(&lock); // same lock as above
   * ... set promise completion...
   * schedule_tid(tid); // wake up the fiber
   * pthread_mutex_unlock(&lock);
   * \endcode
   */
  static void deschedule_self(pthread_mutex_t* lock);

  /**
   *  Schedules a fiber for execution.
   *  If this fiber was previously descheduled by
   *  deschedule_self(), the fiber is scheduled for execution.
   *  Otherwise, nothing happens. Some care must be taken to avoid race
   *  conditions. See the deschedule_self() function for details.
   */
  static void schedule_tid(size_t tid);
};

}

#endif

