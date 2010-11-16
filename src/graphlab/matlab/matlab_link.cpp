#include <mex.h>
#include "gl_emx_graphtypes.hpp"
#include "matlab_link.h"
#include "update_function_generator.hpp"
#include "mx_emx_converters.hpp"

gl_update_function_params* get_params(double handle) {
  // depun the pointer
  // first force cast it back to a uint64_t
  double *handleptr = &handle;
  uint64_t paramsptr = *reinterpret_cast<uint64_t*>(handleptr);
  gl_update_function_params* params;
#ifdef __LP64__
  params = *reinterpret_cast<gl_update_function_params**>(&paramsptr);
  return params;
#else
  uint32_t truncated_paramsptr = paramsptr;
  params = *reinterpret_cast<gl_update_function_params**>(&truncated_paramsptr);
  return params;
#endif
}

void emx_get_edge_data(double handle, uint32_T eid, gl_emx_edgetype *edge) {
  gl_update_function_params *paramsptr = get_params(handle);
  eid--;
  // get the data and copy it out
  const gl_emx_edgetype &e = paramsptr->scope->const_edge_data(eid);
  emxcopy(*edge, e);
}

void emx_get_vertex_data(double handle, uint32_T vid, gl_emx_vertextype *vertex) {
  gl_update_function_params *paramsptr = get_params(handle);
  vid--;
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
  gl_update_function_params *paramsptr = get_params(handle);
  eid--;
  // write the data
  gl_emx_edgetype &e = paramsptr->scope->edge_data(eid);
  emxcopy(e, *edge);
}

void emx_set_vertex_data(double handle, uint32_T vid, gl_emx_vertextype *vertex) {
  gl_update_function_params *paramsptr = get_params(handle);
  vid--;
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


void emx_add_task(double handle, uint32_T vid, const char* fnname, double priority) {
  //std::cout << "add_task: " << vid << " " << fnname << " " << priority << std::endl;
  static bool printed = false;
  gl_update_function_params *paramsptr = get_params(handle);
  vid--;
  if (fnname == NULL) return;
  // figure out the function..
  std::string f(fnname);
  f = "__gl__" + f;
  update_function_map_type::iterator i = update_function_map.find(f);
  if (i != update_function_map.end()) {
    paramsptr->scheduler->add_task(vid, i->second, priority);
  }
  else {
    if (!printed) std::cerr << "Update function " << fnname << " not found." << std::endl;
    printed = true;
  }
}