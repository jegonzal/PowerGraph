#include <mex.h>
#include <graphlab_mex_parse.hpp>
#include <graphlab.hpp>
#include "rtwtypes.h"
#include "updates_types.h"
#include "mx_emx_converters.hpp"
#include "updates_initialize.h"
#include "gl_emx_graphtypes.hpp"


/**
 * graphlab_mex(vertexdata, adj_mat, edgedata, schedule, strict)
 * 
 * vertexdata: cell array of vertex data
 * adj_mat: (sparse) adjacency matrix where adj_mat[i][j] is an edge from vertex
 *          i to vertex j and the data on the edge is edgedata(adjmat[i][j])
 * edgedata: cell array of edge data
 * schedule: array of task structs where each
 *     task struct:
 *       vertex: numeric vertex id
 *       update function: string. Name of the update function
 * strict: numeric. if non-zero, type checks are strictly enforced
 *  Returns 0 on     
 */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // basic data type checks
  // we must output to something
  if (nlhs == 0) {
    mexWarnMsgTxt("No output arguments.");
    plhs[0] = mxCreateDoubleScalar(-1);
    return;
  }
  else if (nlhs > 1) {
    mexWarnMsgTxt("Too many output arguments.");
    return;
  }
  // there must be exactly 4 arguments
  if (nrhs != 5) {
    mexWarnMsgTxt("Erronous function call");
    mexWarnMsgTxt("Usage: graphlab_mex(vertexdata, adj_mat, edgedata, schedule, strict)");
    plhs[0] = mxCreateDoubleScalar(-1);
    return;
  }

  // fill a parameters structure
  mex_parameters param;
  param.vdata = prhs[0];
  param.adjmat = prhs[1];
  param.edata = prhs[2];
  param.sched = prhs[3];
  param.strict = prhs[4];
  // basic type check of the parameters
  if (basic_typecheck(param) == false) {
    plhs[0] = mxCreateDoubleScalar(-1);
    return;
  }
  bool strict = mxGetScalar(param.strict) != 0;
  // make the graph.
  emx_graph graph;
  bool ret = construct_graph(graph, param.vdata, param.adjmat, param.edata);
  if (ret == false) {
    if (strict) {
      mexWarnMsgTxt("Type conversion errors. Strict-mode is set. Terminating.");
      plhs[0] = mxCreateDoubleScalar(-1);
      return;
    }
    else {
      mexWarnMsgTxt("Type conversion errors. Strict-mode is not set. Continuing.");
    }
  }
  updates_initialize();
}
