/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef DC_COMM_SERVICES_HPP
#define DC_COMM_SERVICES_HPP

#include <graphlab/rpc/dc_comm_base.hpp>

namespace graphlab {
  
namespace dc_impl {
/**
 * \ingroup rpc_internal
A set of low level services which operate directly on the comm layer without
going through the RPC controller.
*/
class dc_comm_services {
 private:

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
  
  dc_comm_base *comm;
  
  // set the waiting flag
  void __child_to_parent_barrier_trigger(procid_t source,char releaseval);
  
  void __parent_to_child_barrier_release(char releaseval);
  
 public:
  dc_comm_services(dc_comm_base* comm):comm(comm) { 
    child_barrier[0] = 0; child_barrier[1] = 0;
    barrier_sense = 1;
    barrier_release = 0;
    
    child[0] = (procid_t)(comm->procid() * 2 + 1);
    child[1] = (procid_t)(comm->procid() * 2 + 2); 
    parent =  (procid_t)((comm->procid() - 1) / 2);
  }
 

  void recv(procid_t source, char *c, size_t len);
  void barrier();

};


} // dc_impl
} // graphlab
#endif

