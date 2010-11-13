#ifndef EMLC_LINK_FUNCTIONS_H
#define EMLC_LINK_FUNCTIONS_H

#include "gl_emx_graphtypes.hpp"

void emx_get_edge_data(uint32_T handle, uint32_T eid, gl_emx_edgetype *edge);
void emx_get_vertex_data(uint32_T handle, uint32_T vid, gl_emx_vertextype *vertex);


#endif