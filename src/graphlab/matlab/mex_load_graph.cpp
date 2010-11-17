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
#include "updates_initialize.h"
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
 * vertexdata: cell array of vertex data
 * adj_mat: (sparse) adjacency matrix where adj_mat[i][j] is an edge from vertex
 *          i to vertex j and the data on the edge is edgedata(adjmat[i][j])
 * edgedata: cell array of edge data
 * options: options and schedule struct
 * graphfile: graph output file
 * strict: numeric 0/1 . Strictness of typechecking
 *  Returns new graph  data on exit
 *
 *
 *  optionsstruct:
 * -- scheduler: Scheduler string
 * -- scope: Scope Type
 * -- ncpus: number of cpus to use
 * -- initial_schedule: and array of structs describing the schedule
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
    mexWarnMsgTxt("Usage: graphlab_mex(igraphfile)");
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


