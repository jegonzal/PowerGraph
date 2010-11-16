#ifndef EMLC_LINK_FUNCTIONS_H
#define EMLC_LINK_FUNCTIONS_H

#include "gl_emx_graphtypes.hpp"
/**
 * the matlab get_edge_data wrapper calls this function
 */
void emx_get_edge_data(double handle, uint32_T eid, gl_emx_edgetype *edge);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_get_vertex_data(double handle, uint32_T vid, gl_emx_vertextype *vertex);

/**
 * the matlab set_edge_data wrapper calls this function
 */
void emx_set_edge_data(double handle, uint32_T eid, gl_emx_edgetype *edge);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_set_vertex_data(double handle, uint32_T vid, gl_emx_vertextype *vertex);


#endif