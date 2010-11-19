#include <mex.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <graphlab.hpp>
#include <cstdlib>
#include "graphlab_mex_parse.hpp"
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
 * Parses the options structure
 */
std::vector<parsed_initial_schedule>
parse_emx_schedule(const emxArray_graphlab_initial_schedule &sched) {
  std::vector<parsed_initial_schedule> out;
  size_t numentries = 1;
  for (int i = 0;i < sched.numDimensions; ++i) numentries *= sched.size[i];
  
  for (size_t i = 0;i < numentries; ++i) {

    parsed_initial_schedule parse;
    parse.update_function =
                emxArray_char_T_to_string(sched.data[i].update_function);
    if (parse.update_function.length() == 0) {
      mexWarnMsgTxt("schedule entry with no update function.");
      continue;
    }
    // parse the vertices and priorities
    
    parse.vertices = emxArray_to_vector<uint32_t>(sched.data[i].vertices);
    parse.priorities = emxArray_to_vector<double>(sched.data[i].priorities);
    if (parse.vertices.size() != parse.priorities.size()) {
      mexWarnMsgTxt("#vertices do not match #priorities");
    }
    out.push_back(parse);
  }
  return out;
}



/**
 * mex_save_graph(vertexdata, adj_mat, edgedata, options, graphfile, strict)
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
  if (nlhs != 1) {
    mexWarnMsgTxt("Wrong number of output arguments");
    return;
  }
  if (nrhs != 6) {
    mexWarnMsgTxt("Erronous function call");
    mexWarnMsgTxt("Usage: [.._save_graph](vertexdata, adj_mat, edgedata, options, graphfile, strict)");
    plhs[0] = mxCreateLogicalScalar(false);
    return;
  }

  // fill the parameters structure
  mex_parameters param;
  param.vdata = prhs[0];
  param.adjmat = prhs[1];
  param.edata = prhs[2];
  param.options = prhs[3];
  param.graphfile = prhs[4];
  param.strict = prhs[5];

    // basic type check of the parameters
  if (basic_typecheck(param) == false) {
    mexWarnMsgTxt("Basic typechecks failed.");
    plhs[0] = mxCreateLogicalScalar(false);
    return;
  }
  bool strict = mxGetScalar(param.strict) != 0;

  // read the options information
  graphlab_options_struct optionsstruct;
  memset(&optionsstruct, 0, sizeof(graphlab_options_struct));
  bool ret = mxarray2emx(param.options, optionsstruct);
  if (ret == false) {
    mexWarnMsgTxt("Type conversion errors in Options. Strict-mode is set. Terminating.");
    plhs[0] = mxCreateLogicalScalar(false);
    return;
  }
  // construct the graph
  emx_graph graph;
  ret = construct_graph(graph, param.vdata, param.adjmat, param.edata);
  if (ret == false) {
    if (strict != 0) {
      mexWarnMsgTxt("Type conversion errors in Graph Data. Strict-mode is set. Terminating.");
      cleanup_graph(graph);
      plhs[0] = mxCreateLogicalScalar(false);
      return;
    }
    else {
      mexWarnMsgTxt("Type conversion errors in Graph Data. Strict-mode is not set. Continuing.");
    }
  }
  graph.finalize();
  graph.compute_coloring();

  char* graphfile = mxArrayToString(param.graphfile);
  mexPrintf("Serializing to: %s\n", graphfile);
  std::ofstream fout(graphfile, std::ios_base::binary);
  if (!fout.good()) {
    mexWarnMsgTxt("Unable to open output file! Terminating.");
    cleanup_graph(graph);
    return;
  }
  graphlab::oarchive oarc(fout);
  std::vector<parsed_initial_schedule> schedule = 
                    parse_emx_schedule(*optionsstruct.initial_schedule);
  
  oarc << graph;
  oarc << schedule;
  fout.close();
  cleanup_graph(graph);
  graph.clear();
  plhs[0] = mxCreateLogicalScalar(true);
  return;
}


