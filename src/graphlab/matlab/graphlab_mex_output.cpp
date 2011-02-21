#include <mex.h>
#include "gl_emx_graphtypes.hpp"
#include "mx_emx_converters.hpp"
#include "graphlab_mex_output.hpp"

#include <graphlab/macros_def.hpp>

int output_graph(emx_graph &graph,
                 mxArray* &vdata,
                 mxArray* &adjmat,
                 mxArray* &edata) {

  output_canonical_adj_mat(graph, adjmat);
  // create the vdata cell array
  mwSize dims[2]; dims[0] = 1; dims[1] = graph.num_vertices();
  vdata = mxCreateCellArray(2, dims);

  // output the vdata
  bool conversions_ok = true;
  for (size_t i = 0;i < graph.num_vertices(); ++i) {
    // perform the conversion
    mxArray* mx_vtx = NULL;
    conversions_ok |= emx2mxarray(graph.vertex_data(i), mx_vtx);
    // set the new cell value
    mxSetCell(vdata, i, mx_vtx);
  }

  // create the edata cell array
  dims[0] = 1; dims[1] = graph.num_edges();
  edata = mxCreateCellArray(2, dims);

  // output the edata
  for (size_t i = 0;i < graph.num_edges(); ++i) {
    // perform the conversion
    mxArray* mx_edge = NULL;
    conversions_ok |= emx2mxarray(graph.edge_data(i), mx_edge);
    // set the new value
    mxSetCell(edata, i, mx_edge);
  }
  if (conversions_ok) return 1;
  else return 0;
}



void output_canonical_adj_mat(emx_graph &graph,
                             mxArray* &adjmat) {
  // create a sparse matrix with just the right number of elements
  adjmat = mxCreateSparse(graph.num_vertices(),
                          graph.num_vertices(),
                          graph.num_edges(),
                          mxREAL);

  mwIndex* ir = mxGetIr(adjmat);
  mwIndex* jc = mxGetJc(adjmat);
  double* pr = mxGetPr(adjmat);

  jc[0] = 0;
  size_t idx = 0;
  for (size_t i = 0;i < graph.num_vertices(); ++i) {
    // iterate through all the out edges of vertex i
    foreach(graphlab::edge_id_t eid, graph.out_edge_ids(i)) {
      ir[idx] = graph.target(eid) + 1;
      pr[idx] = eid + 1;
      ++idx;
    }
    jc[i + 1] = idx;
  }
}

