#ifndef DC_SERVICES_HPP
#define DC_SERVICES_HPP
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>

#define DC_SERVICES_BARRIER_BRANCH_FACTOR 128
namespace graphlab {

class dc_services {
 private:
  dc_dist_object<dc_services> rpc;
  
  // Sense reversing barrier
  /// The next value of the barrier. either +1 or -1
  int barrier_sense;
  /// When this flag == the current barrier value. The barrier is complete
  int barrier_release;
  /** when barrier sense is 1, barrier clears when 
   * child_barrier_counter == numchild. When barrier sense is -1, barrier
   * clears when child_barrier_counter == 0;
   */
  atomic<int> child_barrier_counter;
  /// condition variable and mutex protecting the barrier variables
  conditional barrier_cond;
  mutex barrier_mut;
  size_t parent;  /// parent node
  size_t childbase; /// id of my first child
  size_t numchild;  /// number of children
  
  // set the waiting flag
  void __child_to_parent_barrier_trigger(procid_t source);
  
  void __parent_to_child_barrier_release(int releaseval);
  
 public:
  dc_services(distributed_control &dc):rpc(dc, this) { 
    // reset the child barrier values
    child_barrier_counter.value = 0;
    barrier_sense = 1;
    barrier_release = -1;

    // compute my children
    childbase = size_t(rpc.dc().procid()) * DC_SERVICES_BARRIER_BRANCH_FACTOR + 1;
    if (childbase >= rpc.dc().numprocs()) {
      numchild = 0;
    }
    else {
      size_t maxchild = std::min<size_t>(rpc.dc().numprocs(), 
                               childbase + DC_SERVICES_BARRIER_BRANCH_FACTOR);
      numchild = maxchild - childbase;
    }
    
    parent =  (rpc.dc().procid() - 1) / DC_SERVICES_BARRIER_BRANCH_FACTOR   ;
    logger(LOG_INFO, "%d %d %d", childbase, numchild, parent);
  }
 
  /**
   * tree-reduction based sense-reversing barrier with a branching factor of 
   * DC_SERVICES_BARRIER_BRANCH_FACTOR
   */
  void barrier();
  
  
  inline void comm_barrier(procid_t targetmachine) {
    rpc.dc().comm_barrier(targetmachine);
  }
  
  inline void comm_barrier() {
    rpc.dc().comm_barrier();
  }
  
  /**
  This function allows one machine to broadcasts a string of data 
  to all machines.
  
  The originator calls broadcast with data provided in 
  in 'data' and length len. All other machines must call
  broadcast with data = NULL. 
  
  The originator will then return 'data'. All other machines
  will return a new pointer, with the length of the string
  returned in 'len'. The returned pointer must be freed
  by the caller (with the exception of the originator). 
  
  This function is not guaranteed to have barrier-like behavior.
  That is, broadcast could be implemented in a buffered fashion.
  */
  char* broadcast(char* data, size_t &len) { return NULL;};
};

}
#endif
