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

#ifndef COLORING_H
#define COLORING_H

#include <iostream>
#include <graphlab.hpp>

//---------------- TYPES --------------------

/** Numerical representation for color. */
typedef unsigned long color_type;
const unsigned long UNCOLORED = -1;

/**
 * Vertex representation. Each vertex has a color.
 */
struct vertex_data {
  color_type color;
  int saturation;
  vertex_data () : color(UNCOLORED), saturation(0) {};
};

struct edge_data {
  // no edge data required
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;

//--------------- CONSTANTS ---------------
#define DEFAULT_FORMAT      "tsv"
#define NONE                ""

#define OPT_GRAPH_FILE      "graph"
#define OPT_FORMAT          "format"

//--------------- RETURN CODES ------------
#define ERR_INPUT -1
#define SUCCESS 0

#endif
