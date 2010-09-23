#include <graphlab/rpc/dc_services.hpp>

namespace graphlab {
  

void dc_services::__child_to_parent_barrier_trigger(procid_t source) {
  barrier_mut.lock();
  if (source == child[0]) {
    child_barrier[0] = barrier_sense;
  }
  else {
    child_barrier[1] = barrier_sense;
  }
  barrier_cond.signal();
  barrier_mut.unlock();
} 

void dc_services::__parent_to_child_barrier_release(char releaseval) {
  // send the release downwards
  if (child[0] < rpc.dc().numprocs()) {
    rpc.fast_remote_call(child[0],
                        &dc_services::__parent_to_child_barrier_release,
                        releaseval);
  }
  if (child[1] < rpc.dc().numprocs()) {
    rpc.fast_remote_call(child[1],
                        &dc_services::__parent_to_child_barrier_release,
                        releaseval);
  }
  barrier_mut.lock();
  barrier_release = releaseval;    
  barrier_cond.signal();
  barrier_mut.unlock();
}



void dc_services::barrier() {
    // upward message
    char barrier_val = barrier_sense;      
    barrier_mut.lock();

    while(1) {
      if ((child_barrier[0] == barrier_sense || child[0] >= rpc.dc().numprocs()) &&
          (child_barrier[1] == barrier_sense || child[1] >= rpc.dc().numprocs())) {
        // flip the barrier sense
        barrier_sense = ! barrier_sense;
        // call child to parent in parent
        barrier_mut.unlock();
        if (rpc.dc().procid() != 0) {
          rpc.fast_remote_call(parent, 
                              &dc_services::__child_to_parent_barrier_trigger,
                              rpc.dc().procid());
        }
        break;
      }
      barrier_cond.wait(barrier_mut);
    }
    logger(LOG_INFO, "barrier phase 1 complete");
    // I am root. send the barrier releae downwards
    if (rpc.dc().procid() == 0) {
      barrier_release = barrier_val;
      if (child[0] < rpc.dc().numprocs()) {
        rpc.fast_remote_call(child[0],
                            &dc_services::__parent_to_child_barrier_release,
                            barrier_val);
      }
      if (child[1] < rpc.dc().numprocs()) {
        rpc.fast_remote_call(child[1],
                            &dc_services::__parent_to_child_barrier_release,
                            barrier_val);
      } 
    }
    // wait for the downward message releasing the barrier
    barrier_mut.lock();
    while(1) {
      if (barrier_release == barrier_val) break;
      barrier_cond.wait(barrier_mut);
    }
    barrier_mut.unlock();

    logger(LOG_INFO, "barrier phase 2 complete");
  }
  

} // graphlab



