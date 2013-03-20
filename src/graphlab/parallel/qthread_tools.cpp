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
        setenv("QTHREAD_HWPAR", c, 1); 
      }

      qthread_initialize();
      qthread_initialized = true;
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
  retval = 0;
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
  qthread_tools::init();
}

qthread_group::~qthread_group() { 
  join(); 
}

void qthread_group::launch(const boost::function<void (void)> &spawn_function) {
  qthread_thread* thr = new qthread_thread;
  lock.lock();
  threads.push(thr);
  thr->launch(spawn_function);
  lock.unlock();
}

void qthread_group::join() {
  lock.lock();
  while(!threads.empty()) {
    threads.front()->join();
    threads.pop();
  }
  lock.unlock();
}
} // namespace graphlab
