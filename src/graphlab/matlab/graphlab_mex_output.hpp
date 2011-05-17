/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
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

