#ifndef GRAPHLAB_FIBER_CONTROL_HPP
#define GRAPHLAB_FIBER_CONTROL_HPP
#include <cstdlib>
#include <boost/context/all.hpp>
#include <boost/function.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
namespace graphlab {

/**
 * The master controller for the user mode threading system
 */
class fiber_control {
 public:
  struct fiber {
    simple_spinlock lock;
    fiber_control* parent;
    boost::context::fcontext_t* context;
    void* stack;
    size_t id;
    int affinity;
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
    bool priority;  // flag. If set, rescheduling this fiber
                    // will cause it to be placed at the head of the queue
  };


 private:
  size_t nworkers;
  atomic<size_t> fiber_id_counter;
  atomic<size_t> fibers_active;
  mutex join_lock;
  conditional join_cond;

  bool stop_workers;

  // The scheduler is a simple queue. One for each worker
  struct thread_schedule {
    mutex active_lock;
    conditional active_cond;
    // used as a sentinel for the actual queues
    // first element is always head,
    // tail points to the last element
    fiber active_head;
    fiber* active_tail;
    size_t nactive;
    volatile bool waiting;
  };
  std::vector<thread_schedule> schedule;

  atomic<size_t> nactive;
  thread_group workers;


  // locks must be acquired outside the call
  void active_queue_insert_head(size_t workerid, fiber* value);
  void active_queue_insert_tail(size_t workerid, fiber* value);
  fiber* active_queue_remove(size_t workerid);

  // a thread local storage for the worker to point to a fiber
  static bool tls_created;
  struct tls {
    fiber_control* parent;
    fiber* prev_fiber; // the fiber we context switch from
    fiber* cur_fiber; // the fiber we are context switching to
    fiber* garbage; // A fiber to delete after the context switch
    size_t workerid;
    boost::context::fcontext_t base_context;
  };

  static pthread_key_t tlskey; // points to the tls structure above
  static void tls_deleter(void* tls);

  /// internal function to create the TLS for the worker threads
  static void create_tls_ptr();
  /// internal function to read the TLS for the worker threads
  static tls* get_tls_ptr();
  /// Returns the current fiber scheduled on this worker thread
  static fiber* get_active_fiber();

  /// The function that each worker thread starts off running
  void worker_init(size_t workerid);

  void reschedule_fiber(size_t workerid, fiber* pfib);
  void yield_to(fiber* next_fib);
  static void trampoline(intptr_t _args);

  size_t load_balanced_worker_choice(size_t seed);

  void (*flsdeleter)(void*);
  
 public:

  /// Private constructor
  fiber_control(size_t nworkers, size_t affinity_base);

  ~fiber_control();

  /** the basic launch function
   * Returns a fiber ID. IDs are not sequential.
   * \note The ID is really a pointer to a fiber_control::fiber object.
   */
  size_t launch(boost::function<void (void)> fn, 
                size_t stacksize = 8192, 
                int worker_affinity = -1);

  inline size_t num_active() const {
    return nactive;
  }
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
   * Implodes dramatically if called from outside a fiber.
   */
  static void exit();

  /**
   * Yields to another fiber.
   * Note that this function will only work within a fiber.
   * If called from outside a fiber, returns immediately.
   */
  static void yield();


  /// True if the singleton instance was created
  static bool instance_created; 
  static size_t instance_construct_params_nworkers; 
  static size_t instance_construct_params_affinity_base;

  /**
   * Sets the fiber control construction parameters.
   * Fails with an assertion failure if the instance has already been created.
   * Must be called prior to any other calls to get_instance()
   * \param nworkers Number of worker threads to spawn. If set to 0,
   *                 the number of workers will be automatically determined
   *                 based on the number of cores the system has.
   * \param affinity_base First worker will have CPU affinity equal to 
   *                      affinity_base. Second will be affinity_base + 1, etc.
   *                      Defaults to 0.
   */
  static void instance_set_parameters(size_t nworkers,
                                      size_t affinity_base);

  /**
   * Gets a reference to the main fiber control singleton
   */
  static fiber_control& get_instance();

  /**
   * Returns the current fiber handle.
   * Note that fiber handles are not sequential, and are really a
   * pointer to an internal datastructure.
   * If called from within a fiber, returns a non-zero value.
   * If called out outsize a fiber, returns 0.
   */
  static size_t get_tid();


  /**
   * Returns true if the calling thread is in a fiber, false otherwise.
   */
  static bool in_fiber();

  /**
   * Returns the worker managing the current fiber.
   * Worker IDs are sequential.
   * If called from outside a fiber, returns (size_t)(-1)
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
   *  This thread by default will be stuck at the head of queue
   *  and will wake up quickly.
   *
   *  \param priority If true, thread will be placed at the head
   *  of the scheduler. If false, it will be placed at the tail
   *  of the scheduler
   */
  static void schedule_tid(size_t tid, bool priority = true);
};

}

#endif

