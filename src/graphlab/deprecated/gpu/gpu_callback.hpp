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


#ifndef GRAPHLAB_GPU_CALLBACK_HPP
#define GRAPHLAB_GPU_CALLBACK_HPP

#include <vector>


#include <graphlab/graph/graph.hpp>


namespace graphlab {

  class gpu_callback {
  public:

    /**
     * Adds a task to execute the update function on the vertex with
     * the given priority.
     */
    void add_task(vertex_id_t vertex, double priority) { }
                         
  };

}; //end graphlab namespace

#endif

