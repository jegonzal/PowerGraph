#include <graphlab/parallel/qthread_tools.hpp>
#include <cstdlib>
#include <cstdio>

namespace graphlab {

namespace qthread_tools {
  void init(int numworkers, int stacksize) {
    static bool qthread_initialized = false;
    if (!qthread_initialized) {
      if (stacksize > 0) {
        // we need to set the environment variable to force the stacksize
        char c[32];
        sprintf(c, "%d", stacksize);
        setenv("QTHREAD_STACK_SIZE", c, 1);
      }
      if (numworkers > 0) {
        // we need to set the environment variable to force the number of workers
        char c[32];
        sprintf(c, "%d", numworkers);
        //setenv("QTHREAD_HWPAR", c, 1);
        setenv("QTHREAD_NUM_WORKERS_PER_SHEPHERD", c, 1);
      }
      // force the number of shepherds to 1 so all threads can migrate
      {
        char c[32];
        sprintf(c, "%d", 1);
        setenv("QTHREAD_NUM_SHEPHERDS", c, 1);
      }
      qthread_initialize();
      qthread_initialized = true;
      std::cout << "Number of workers: " << qthread_num_workers() << "\n";
      std::cout << "Number of shepherds: " << qthread_num_shepherds() << "\n";
    }
  }


  void finalize() {
    static bool qthread_finalized = false;
    if (!qthread_finalized) {
      qthread_finalize();
      qthread_finalized = true;
    }
  }
} // qthread_tools




//----------- qthread_thread ---------------
qthread_thread::qthread_thread() {
  qthread_tools::init();
  retval = 0;
}

qthread_thread::~qthread_thread() {
  join();
}

aligned_t qthread_thread::invoke(void *_args) {
  qthread_thread::invoke_args* args =
      reinterpret_cast<qthread_thread::invoke_args*>(_args);
  args->spawn_routine();
  delete args;
  return 0;
}


void qthread_thread::launch(const boost::function<void (void)> &spawn_routine) {
  qthread_fork(qthread_thread::invoke,
               static_cast<void*>(new invoke_args(spawn_routine)),
               &retval);
}

void qthread_thread::join() {
  qthread_readFF(NULL, &retval);
}



void qthread_thread::yield() {
  qthread_yield();
}

//--------------- qthread_group ------------------

qthread_group::qthread_group() {
}

qthread_group::~qthread_group() {
  join();
}

void qthread_group::launch(const boost::function<void (void)> &spawn_function) {
  qthread_thread* thr = new qthread_thread;
  thr->launch(spawn_function);
  lock.lock();
  threads.push(thr);
  lock.unlock();
}

void qthread_group::join() {
  lock.lock();
  while(!threads.empty()) {
    qthread_thread* thread_ptr = threads.front();
    lock.unlock();
    thread_ptr->join();
    delete thread_ptr;
    lock.lock();
    threads.pop();
  }
  lock.unlock();
}
} // namespace graphlab
