// unfortunately this file has to be .h for emlc to recognize it
#ifndef EMLC_LINK_FUNCTIONS_H
#define EMLC_LINK_FUNCTIONS_H

#include "gl_emx_graphtypes.hpp"
/**
 * the matlab get_edge_data wrapper calls this function
 */
void emx_get_edge_data(HANDLE_TYPE handle, uint32_T eid, gl_emx_edgetype *edge);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_get_vertex_data(HANDLE_TYPE handle, uint32_T vid, gl_emx_vertextype *vertex);

/**
 * the matlab set_edge_data wrapper calls this function
 */
void emx_set_edge_data(HANDLE_TYPE handle, uint32_T eid, gl_emx_edgetype *edge);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_set_vertex_data(HANDLE_TYPE handle, uint32_T vid, gl_emx_vertextype *vertex);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_set_vertex_data(HANDLE_TYPE handle, uint32_T vid, gl_emx_vertextype *vertex);

/**
 * The matlab add_task wrapper calls this function
 */
void emx_add_task(HANDLE_TYPE handle, uint32_T vid, const char* fnname, double priority);
#endif