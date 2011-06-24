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


#ifndef GRAPHLAB_MEX_TYPECHECK_HPP
#define GRAPHLAB_MEX_TYPECHECK_HPP
#include "mex_save_graph.hpp"
#include "gl_emx_graphtypes.hpp"

/**
 * Performs a check over all the arguments to ensure that the
 * basic input formats are satisfied.
 * Outputs format failures as Matlab warnings.
 * Returns false on failure, true on success.
 */
bool basic_typecheck(mex_parameters &param);

/**
 * Constructs a graph from the vdata, edata, adjmat parameters
 * Return 1 on success.
 * Returns 0 on non-fatal conversion errors/warnings
 * Return -1 on failures
 */
int construct_graph(emx_graph &graph, const mxArray* vdata,
                    const mxArray* adjmat, const mxArray* edata);

/**
 * Tries to insert the edge from src->dest using the data in edata[edataidx]
 * Return 1 on success.
 * Returns 0 on non-fatal conversion errors/warnings
 * Return -1 on failures
 */
int try_add_edge(emx_graph &graph, const mxArray* edata,
                 size_t src, size_t dest, size_t edataidx);

/**
 * Inserts all edges described in adjmat into the graph.
 * adjmat must be a full matrix.
 */
int add_edges_from_full_adj(emx_graph &graph,
                            const mxArray* adjmat,
                            const mxArray* edata);

/**
 * Inserts all edges described in adjmat into the graph.
 * adjmat must be a sparse matrix.
 */
int add_edges_from_sparse_adj(emx_graph &graph,
                              const mxArray* adjmat,
                              const mxArray* edata);
#endif

