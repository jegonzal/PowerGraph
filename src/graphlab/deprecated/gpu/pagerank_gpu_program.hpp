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


#ifndef GRAPHLAB_PAGERANK_GPU_PROGRAM
#define GRAPHLAB_PAGERANK_GPU_PGOGRAM


namespace graphlab {

  class pagerank_gpu_program {
  public:
    // Struct for vertex data
    struct vertex_data_type {
      float value;
      float self_weight;
    };

    // struct for edge data
    struct edge_data_type {
      float weight;
      float old_source_value;
    };

    typedef gpu_scope<pagerank_gpu_program> scope_type;
    typedef gpu_callback                    callback_type;

    __device__ static void 
    update_function(scope_type& scope,
                    callback_type& callback) {
                    
      // Get the data associated with the vertex
      vertex_data& vdata = scope.vertex_data();
  
      // Sum the incoming weights; start by adding the 
      // contribution from a self-link.
      float sum = vdata.value * vdata.self_weight;
  
      foreach(graphlab::edge_id_t eid, scope.in_edge_ids()) {
        // Get the neighobr vertex value
        const vertex_data& neighbor_vdata =
          scope.const_neighbor_vertex_data(scope.source(eid));
        double neighbor_value = neighbor_vdata.value;
    
        // Get the edge data for the neighbor
        edge_data& edata = scope.edge_data(eid);
        // Compute the contribution of the neighbor
        double contribution = edata.weight * neighbor_value;
    
        // Add the contribution to the sum
        sum += contribution;
    
        // Remember this value as last read from the neighbor
        edata.old_source_value = neighbor_value;
      }

      // compute the jumpweight
      sum = (1-damping_factor)/scope.num_vertices() + damping_factor*sum;
      vdata.value = sum;
  
      // Schedule the neighbors as needed
      foreach(graphlab::edge_id_t eid, scope.out_edge_ids()) {
        edge_data& outedgedata = scope.edge_data(eid);
    
        // Compute edge-specific residual by comparing the new value of this
        // vertex to the previous value seen by the neighbor vertex.
        double residual =
          outedgedata.weight *
          std::fabs(outedgedata.old_source_value - vdata.value);
        // If the neighbor changed sufficiently add to scheduler.
        if(residual > termination_bound) {
          gl_types::update_task task(scope.target(eid), pagerank_update);
          scheduler.add_task(task, residual);
        }
      }



    } // end of update_function
    
    

  }; // end of gpu pagerank program



};



#endif

