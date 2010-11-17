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