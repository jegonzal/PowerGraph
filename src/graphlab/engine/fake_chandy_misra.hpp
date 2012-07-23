/*  
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


#ifndef GRAPHLAB_FAKE_CHANDY_MISRA_HPP
#define GRAPHLAB_FAKE_CHANDY_MISRA_HPP
#include <vector>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/engine/chandy_misra_interface.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

/**
  * 
  * \internal
  */
template <typename GraphType>
class fake_chandy_misra: public chandy_misra_interface<GraphType> {
 public:
  typedef typename GraphType::local_vertex_type local_vertex_type;
  typedef typename GraphType::local_edge_type local_edge_type;
  
  typedef typename GraphType::vertex_id_type vertex_id_type;
  typedef typename GraphType::lvid_type lvid_type;

  typedef fake_chandy_misra<GraphType> dcm_type;
 private:
  GraphType &graph;
  dc_dist_object<dcm_type> rmi;
  std::vector<atomic<procid_t> > counter;

  boost::function<void(lvid_type)> callback;
  boost::function<void(lvid_type)> hors_doeuvre_callback;

  void rpc_decrement_counter(vertex_id_type gvid) {
    decrement_counter(graph.vertex(gvid).local_id());
  }
  void decrement_counter(lvid_type lvid) {
    if (counter[lvid].dec() == 0) {
      callback(lvid);
    }
  }
 public:
  inline fake_chandy_misra(distributed_control &dc,
                                  GraphType &graph,
                                  boost::function<void(lvid_type)> callback,
                                  boost::function<void(lvid_type)> hors_doeuvre_callback = NULL
                                  ):
                          rmi(dc, this),
                          graph(graph),
                          callback(callback),
                          hors_doeuvre_callback(hors_doeuvre_callback){
    counter.resize(graph.num_local_vertices());
    for (lvid_type i = 0;i < graph.num_local_vertices(); ++i) {
      local_vertex_type lvertex(graph.l_vertex(i));
      if (lvertex.owner() == rmi.procid()) {
        counter[i] = lvertex.num_mirrors() + 1;
      }
      else {
        counter[i] = 0;
      }
    }
    rmi.barrier();
  }

  size_t num_clean_forks() const {
    return 0; 
  }

  void make_philosopher_hungry_per_replica(lvid_type p_id) {
    hors_doeuvre_callback(p_id);
    // signal master to complete
    // get master
    local_vertex_type lvertex(graph.l_vertex(p_id));
    procid_t owner = lvertex.owner();
    if (owner == rmi.procid()) decrement_counter(p_id);
    else rmi.remote_call(owner, &dcm_type::rpc_decrement_counter, lvertex.global_id());
  }
  
  
  void philosopher_stops_eating_per_replica(lvid_type p_id) {
   if (graph.l_vertex(p_id).owner() == rmi.procid()) {
     counter[p_id] = graph.l_vertex(p_id).num_mirrors() + 1;
    }
  }
};

}

#include <graphlab/macros_undef.hpp>

#endif
