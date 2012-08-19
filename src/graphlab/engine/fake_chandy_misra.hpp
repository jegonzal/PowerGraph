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
  boost::function<void(lvid_type)> callback;
  boost::function<void(lvid_type)> hors_doeuvre_callback;

 public:
  inline fake_chandy_misra(distributed_control &dc,
                                  GraphType &graph,
                                  boost::function<void(lvid_type)> callback,
                                  boost::function<void(lvid_type)> hors_doeuvre_callback = NULL
                                  ):
                          callback(callback),
                          hors_doeuvre_callback(hors_doeuvre_callback){
  }

  size_t num_clean_forks() const {
    return 0; 
  }

  void make_philosopher_hungry_per_replica(lvid_type p_id) {
  }
  
  
  void philosopher_stops_eating_per_replica(lvid_type p_id) {
  }
};

}

#include <graphlab/macros_undef.hpp>

#endif
