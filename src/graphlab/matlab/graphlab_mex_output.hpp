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


#ifndef GRAPHLAB_MEX_OUTPUT_HPP
#define GRAPHLAB_MEX_OUTPUT_HPP
#include "gl_emx_graphtypes.hpp"

/**
 * Outputs the vdata and edata of the graph. vdata and edata should be empty
 * Return 1 on success.
 * Returns 0 on non-fatal conversion errors/warnings
 * Return -1 on failures
 */
int output_graph(emx_graph &graph, 
                 mxArray* &vdata,
                 mxArray* &adjmat,
                 mxArray* &edata);


/**
 * Outputs the graph structure of the graph as a sparse matrix ,
 * where each entry is the edge id adjmat must be empty.
 */
void output_canonical_adj_mat(emx_graph &graph,
                             mxArray* &adjmat);
#endif

