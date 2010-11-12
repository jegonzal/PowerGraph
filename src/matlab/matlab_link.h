#ifndef EMLC_LINK_FUNCTIONS_H
#define EMLC_LINK_FUNCTIONS_H

#include <tmwtypes.h>

void get_in_edges_impl(int32_T handle, int32_T vertex, emxArray_int32_T *in_edges);

void get_out_edges_impl(int32_T handle, int32_T vertex, emxArray_int32_T *out_edges);

int32_T get_src_vertex_impl(int32_T handle, int32_T edge);
int32_T get_dest_vertex_impl(int32_T handle, int32_T edge);


#endif