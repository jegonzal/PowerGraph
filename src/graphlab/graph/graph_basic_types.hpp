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


#ifndef GRAPHLAB_GRAPH_BASIC_TYPES
#define GRAPHLAB_GRAPH_BASIC_TYPES

#include <stdint.h>

namespace graphlab {
  /// Identifier type of a vertex which is globally consistent.
  typedef uint32_t vertex_id_type;
  
  /// Identifier type of a vertex which is only locally consistent.
  typedef uint32_t lvid_type;
  
  /// Identifier type of an edge which is only locally consistent.
  typedef uint32_t edge_id_type;
  
  /**
   * The set of edges that are operated on during gather and scatter
   * operations.
   */
  enum edge_dir_type {NO_EDGES = 0, IN_EDGES = 1, OUT_EDGES = 2 , ALL_EDGES = 3};
} // end of namespace graphlab

#endif
