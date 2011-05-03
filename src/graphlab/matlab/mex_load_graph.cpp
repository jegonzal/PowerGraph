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

#include <mex.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <graphlab.hpp>
#include <cstdlib>
#include "graphlab_mex_output.hpp"
#include "rtwtypes.h"
#include "updates_types.h"
#include "mx_emx_converters.hpp"
#include "gl_emx_graphtypes.hpp"

/**
 * Frees the graph datastructures.
 */
void cleanup_graph(emx_graph &graph) {
  // free the graph
  for (size_t i = 0;i < graph.num_vertices(); ++i) {
    freeemx(graph.vertex_data(i));
  }
  for (size_t i = 0;i < graph.num_edges(); ++i) {
    freeemx(graph.edge_data(i));
  }
  graph.clear();
}


/**
 * [vertexdata, adj_mat, edgedata] mex_load_graph(graphfile)
*
* graphfile: graph output file
 * strict: numeric 0/1 . Strictness of typechecking
 *  Returns new graph  data on exit
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // basic data type checks
  // there must be exactly 6 arguments
  if (nlhs != 3) {
    mexWarnMsgTxt("Not the right number of output arguments");
    return;
  }
  if (nrhs != 1) {
    mexWarnMsgTxt("Erronous function call");
    mexWarnMsgTxt("Usage: [.._load_graph](igraphfile)");
    return;
  }

  // fill the parameters structure
  const mxArray* inputgraphfile = prhs[0];
  if (inputgraphfile == NULL || !mxIsChar(inputgraphfile)) {
    mexWarnMsgTxt("Input should be a file name.");
    return;
  }
  // construct the graph
  emx_graph graph;
  char* graphfile = mxArrayToString(inputgraphfile);
  
  mexPrintf("Deserializing from: %s\n", graphfile);
  std::ifstream fin(graphfile, std::ios_base::binary);
  if (!fin.good()) {
    mexWarnMsgTxt("Unable to open input file! Terminating.");
    return;
  }
  graphlab::iarchive iarc(fin);
  iarc >> graph;
  fin.close();
  output_graph(graph, plhs[0], plhs[1], plhs[2]);

  cleanup_graph(graph);
  graph.clear();
  return;
}


