#include <graphlab/rpc/dc_services.hpp>
#include <graphlab/logger/assertions.hpp>
namespace graphlab {
  

void dc_services::__child_to_parent_barrier_trigger(procid_t source) {
  barrier_mut.lock();
  // assert childbase <= source <= childbase + DC_SERVICES_BARRIER_BRANCH_FACTOR
  ASSERT_GE(source, childbase);
  ASSERT_LT(source, childbase + DC_SERVICES_BARRIER_BRANCH_FACTOR);
  child_barrier_counter.inc(barrier_sense);
  barrier_cond.signal();
  barrier_mut.unlock();
} 

void dc_services::__parent_to_child_barrier_release(int releaseval) {
  // send the release downwards
  // get my largest child
  for (size_t i = 0;i < numchild; ++i) {
    rpc.fast_remote_call(childbase + i,
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
  // wait for all children to be done
  while(1) {
    if ((barrier_sense == -1 && child_barrier_counter.value == 0) || 
        (barrier_sense == 1 && child_barrier_counter.value == (int)(numchild))) {
      // flip the barrier sense
      barrier_sense = -barrier_sense;
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
  // I am root. send the barrier release downwards
  if (rpc.dc().procid() == 0) {
    barrier_release = barrier_val;

    for (size_t i = 0;i < numchild; ++i) {
      rpc.fast_remote_call(childbase + i,
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



