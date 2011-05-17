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
#include "gl_emx_graphtypes.hpp"
#include "mx_emx_converters.hpp"
#include "graphlab_mex_parse.hpp"


bool basic_typecheck(mex_parameters &param) {
  // vdata is a cell array of vertex data
  if (param.vdata == NULL || !mxIsCell(param.vdata)) {
    mexWarnMsgTxt("First parameter should be a cell array of vertex data.");
    return false;
  }

  // adjmat is a numeric matrix
  if (param.adjmat == NULL || !mxIsDouble(param.adjmat)) {
    mexWarnMsgTxt("Second parameter should be a (sparse) numeric (double) adjacency matrix.");
  }

  
  // edata is a cell array of edge data
  if (param.edata == NULL || !mxIsCell(param.edata)) {
    mexWarnMsgTxt("Third parameter should be a cell array of edge data.");
    return false;
  }

  // options is a struct array
  if (param.options == NULL || !mxIsStruct(param.options)) {
    mexWarnMsgTxt("Fourth parameter should be a struct array describing options.");
    return false;
  }
  // graphfile is a filename
  if (param.graphfile == NULL || !mxIsChar(param.graphfile)) {
    mexWarnMsgTxt("Fifth parameter should be a file name.");
    return false;
  }

  // strict field should be a string
  if (param.strict == NULL || !mxIsNumeric(param.strict)) {
    mexWarnMsgTxt("strict should be a numeric scalar.");
    return false;
  }
  return true;
}


int construct_graph(emx_graph &graph,
                    const mxArray* vdata,
                    const mxArray* adjmat,
                    const mxArray* edata) {
  size_t nvertices = mxGetM(vdata) * mxGetN(vdata);
  mexPrintf("Vertices: %d\n", nvertices);
  // add vertices
  bool conversions_ok = true;
  for (size_t i = 0;i < nvertices; ++i) {
    // get the cell entry
    const mxArray* mx_vtx = mxGetCell(vdata, i);
    gl_emx_vertextype emx_vtx;
    if (mx_vtx == NULL) {
      clearemx(emx_vtx);
      // mx_vtx is empty. insert an empty entry. and
      // denote as conversion failure
      graph.add_vertex(emx_vtx);
      conversions_ok = false;
    }
    else {
      // perform the conversion
      conversions_ok |= mxarray2emx(mx_vtx, emx_vtx);
      // insert the vertex
      graph.add_vertex(emx_vtx);
    }
  }

  // add edges
  int ret;
  if (mxIsSparse(adjmat)) {
    ret = add_edges_from_sparse_adj(graph, adjmat, edata);
  }
  else {
    ret = add_edges_from_full_adj(graph, adjmat, edata);
  }
  if (ret < 0) return -1;

  conversions_ok |= (ret == 1);
  
  if (conversions_ok) return 1;
  else return 0;
  
}


int try_add_edge(emx_graph &graph, const mxArray* edata,
                 size_t src, size_t dest, size_t edataidx) {
  if (src == dest) {
    char errmsg[256];
    sprintf(errmsg, "Error: cannot add self edge %d -> %d", int(edataidx + 1), int(edataidx + 1));
    mexWarnMsgTxt(errmsg);
    return -1;
  }
  // get the edge data
  const mxArray* mx_edge = mxGetCell(edata, edataidx);

    // see whether the edge data exists
  if (mx_edge == NULL) {
    // edge data out of range!
    char errmsg[256];
    sprintf(errmsg, "Error: edge data entry %d does not exist", int(edataidx + 1));
    mexWarnMsgTxt(errmsg);
    return -1;
  }
  // convert
  gl_emx_edgetype emx_edge;
  bool conversion_ok = mxarray2emx(mx_edge, emx_edge);
  
  // insert the edge
  // but first make sure that source and dest are in range.
  if (src < graph.num_vertices() && dest < graph.num_vertices()) {
    graph.add_edge(src, dest, emx_edge);
  }
  else {
    char errmsg[256];
    sprintf(errmsg, "Error: Adding edge (%d,%d) when there are only %d vertices",
                    int(src + 1), int(dest + 1), int(graph.num_vertices()));
    mexWarnMsgTxt(errmsg);
    return -1;
  }

  if (conversion_ok) return 1;
  else return 0;
}

int add_edges_from_full_adj(emx_graph &graph,
                            const mxArray* adjmat,
                            const mxArray* edata) {
  bool conversions_ok = true;
  size_t m = mxGetM(adjmat);
  size_t n = mxGetN(adjmat);
  double* ptr = mxGetPr(adjmat);
  // no edges
  if (ptr == NULL) {
    return 1;
  }
  for (size_t i_n = 0; i_n < n; ++i_n) {
    for (size_t i_m = 0; i_m < m; ++i_m) {
      size_t edataidx = (*ptr);
      if (edataidx > 0) {
        edataidx--; // 0 base
        int ret = try_add_edge(graph, edata, i_n, i_m, edataidx);
        if (ret < 0) return -1;
        conversions_ok |= (ret ==1);
      }
      ptr++;
    }
  }

  mexPrintf("Edges: %d\n", graph.num_edges());
  if (conversions_ok) return 1;
  else return 0;
}

int add_edges_from_sparse_adj(emx_graph &graph,
                              const mxArray* adjmat,
                              const mxArray* edata) {
  bool conversions_ok = true;
  // this part is hard to explain. Read the documentation for mxGetIr,
  // and mxGetJc
  mwIndex* ir = mxGetIr(adjmat);
  mwIndex* jc = mxGetJc(adjmat);
  double* pr = mxGetPr(adjmat);
  int nedges = mxGetNzmax(adjmat);


  mexPrintf("Edges: %d\n", nedges);
  
  int numcols = mxGetN(adjmat);
  // look for non-zero elements
  int i = 0;
  int j = 0;
  int numinthiscol = jc[1] - jc[0];
  for (int idx = 0; idx < nedges; ++idx) {
    while (numinthiscol <= 0 && j < numcols - 1) {
      j++;
      numinthiscol = jc[j + 1] - jc[j];
    }
    if (numinthiscol == 0) continue;
    i = ir[idx];
    --numinthiscol;
    if (pr[idx] > 0) {
      size_t edataidx = pr[idx] - 1;  // arrays now are 0-based so -1
      int ret = try_add_edge(graph, edata, i, j, edataidx);
      if (ret < 0) return -1;
      conversions_ok |= (ret ==1);
    }
  }

  if (conversions_ok) return 1;
  else return 0;
}

