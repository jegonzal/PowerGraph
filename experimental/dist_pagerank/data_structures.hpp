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


#ifndef PAGERANK_DATA_STRUCTURES
#define PAGERANK_DATA_STRUCTURES

#include <cstring>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <sstream>


#include <graphlab.hpp>
#include <distributed_graphlab.hpp>





/**
 * Stores the value and the self weight
 */
struct vertex_data {
  float value;
  // we save weight of vertex's self-edge in the vertex data
  vertex_data(float value = 1) : value(value) { }
}; // End of vertex data
SERIALIZABLE_POD(vertex_data);

typedef char edge_data;




typedef graphlab::core<vertex_data, edge_data> core_type;
typedef core_type::types::graph local_graph_type;

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::distributed_types<graph_type> gl_types;



#endif
