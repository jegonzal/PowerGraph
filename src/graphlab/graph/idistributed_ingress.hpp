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

#ifndef GRAPHLAB_IDISTRIBUTED_INGRESS_HPP
#define GRAPHLAB_IDISTRIBUTED_INGRESS_HPP

#include <vector>
#include <graphlab/graph/graph_basic_types.hpp>


namespace graphlab {

  template<typename VertexData, typename EdgeData>
  class idistributed_ingress {
  public:     
    virtual void add_edge(vertex_id_type source, vertex_id_type target,
                          const EdgeData& edata) = 0;
    virtual void add_vertex(vertex_id_type vid, const VertexData& vdata) = 0;
    virtual void finalize() = 0;
  }; // end of idstributed_ingress

}; // end of namespace graphlab



#endif
