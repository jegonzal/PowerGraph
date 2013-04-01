#ifndef GRAPHLAB_FIBER_HPP
#define GRAPHLAB_FIBER_HPP
#include <cstdlib>
#include <boost/context/all.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
namespace graphlab {

/// a very simple user-mode threading
class fiber_group {
 public:
  struct fiber {
    boost::context::fcontext_t* context;
    fiber_group* parent;
    void* stack;
    size_t id;
    void* tls;
    fiber* next;
    intptr_t initial_trampoline_args;
    bool terminate;
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
    boost::context::fcontext_t base_context;
  };
  static pthread_key_t tlskey;
  static void tls_deleter(void* tls);

  static void create_tls_ptr();
  static tls* get_tls_ptr();
  static fiber* get_active_fiber();
  static void set_active_fiber(fiber*);

  void worker_init();

  void reschedule_fiber(fiber* pfib);
  void yield_to(fiber* next_fib);
  static void trampoline(intptr_t _args);

 public:
  /// creates a group of fibers using a certain number of worker threads.
  fiber_group(size_t nworkers, size_t stacksize);

  ~fiber_group();
  /** the basic launch function
   * More advanced ones (for instance using boost::function)
   * can be built on top of this easily.
   * Returns a fiber ID. IDs are not sequential.
   * \note The ID is really a pointer to a fiber_group::fiber object.
   */
  size_t launch(void fn(void*), void* param);

  void join();
  static void exit();
  static void yield();
};

}

#endif

