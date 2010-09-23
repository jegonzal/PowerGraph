#ifndef DC_SERVICES_HPP
#define DC_SERVICES_HPP
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_services_base.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>

namespace graphlab {
class dc_services:public dc_services_base {
 private:
  dc_dist_object<dc_services> rpc;
  //DC_DIST_OBJECT_ACCESS(dc_services);
  // barrier flags. 
  /// The next value of the barrier 
  char barrier_sense;
  /// When this flag == the current barrier value. The barrier is complete
  char barrier_release;
  /// Set to the current value of the barrier when the child is done
  char child_barrier[2];
  /// condition variable and mutex protecting the barrier variables
  conditional barrier_cond;
  mutex barrier_mut;
  procid_t parent;  /// parent node
  procid_t child[2]; /// children nodes
  
  
  // set the waiting flag
  void __child_to_parent_barrier_trigger(procid_t source);
  
  void __parent_to_child_barrier_release(char releaseval);
  
 public:
  dc_services(distributed_control &dc):rpc(dc, this) { 
    child_barrier[0] = 0; child_barrier[1] = 0;
    barrier_sense = 1;
    barrier_release = 0;
    
    child[0] = rpc.dc().procid() * 2 + 1;
    child[1] = rpc.dc().procid() * 2 + 2; 
    parent =  (rpc.dc().procid() - 1) / 2;
  }
 

  void barrier();
  
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
