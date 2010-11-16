#include <mex.h>
#include "gl_emx_graphtypes.hpp"
#include "matlab_link.h"
#include "update_function_generator.hpp"
#include "mx_emx_converters.hpp"

void emx_get_edge_data(double handle, uint32_T eid, gl_emx_edgetype *edge) {
  // turn the handle back into a pointer
  double *handleptr = &handle;
  gl_update_function_params *paramsptr = reinterpret_cast<gl_update_function_params*>(handleptr);
  // get the data and copy it out
  const gl_emx_edgetype &e = paramsptr->scope->const_edge_data(eid);
  emxcopy(*edge, e);
}

void emx_get_vertex_data(double handle, uint32_T vid, gl_emx_vertextype *vertex) {
  // turn the handle back into a pointer
  double *handleptr = &handle;
  gl_update_function_params *paramsptr = reinterpret_cast<gl_update_function_params*>(handleptr);
  // get the data and copy it out
  // if vid is current vertex, we use vertex_data. Otherwise we use neighbor_vertex_data
  if (vid == paramsptr->scope->vertex()) {
    const gl_emx_vertextype &v = paramsptr->scope->vertex_data();
    emxcopy(*vertex, v);
  }
  else {
    const gl_emx_vertextype &v = paramsptr->scope->const_neighbor_vertex_data(vid);
    emxcopy(*vertex, v);
  }
}

void emx_set_edge_data(double handle, uint32_T eid, gl_emx_edgetype *edge) {
  // turn the handle back into a pointer
  double *handleptr = &handle;
  gl_update_function_params *paramsptr = reinterpret_cast<gl_update_function_params*>(handleptr);
  // write the data
  gl_emx_edgetype &e = paramsptr->scope->edge_data(eid);
  emxcopy(e, *edge);
}

void emx_set_vertex_data(double handle, uint32_T vid, gl_emx_vertextype *vertex) {
  // turn the handle back into a pointer
  double *handleptr = &handle;
  gl_update_function_params *paramsptr = reinterpret_cast<gl_update_function_params*>(handleptr);
  // write the data
  // if vid is current vertex, we use vertex_data. Otherwise we use neighbor_vertex_data
  if (vid == paramsptr->scope->vertex()) {
    gl_emx_vertextype &v = paramsptr->scope->vertex_data();
    emxcopy(v, *vertex);
  }
  else {
    gl_emx_vertextype &v = paramsptr->scope->neighbor_vertex_data(vid);
    emxcopy(v, *vertex);
  }
}
