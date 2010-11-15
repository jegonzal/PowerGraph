#include <mex.h>
#include <graphlab.hpp>
#include "rtwtypes.h"
#include "updates_types.h"
#include "mx_emx_converters.hpp"
#include "updates_initialize.h"
#include "updates.h"


/**
 * graphlab_mex(vertexdata, adj_mat, edgedata, schedule)
 * vertexdata: cell array of vertex data
 * adj_mat: sparse adjacency matrix where adj_mat[i][j] is an edge from vertex
 *          i to vertex j and the data on the edge is edgedata(adjmat[i][j])
 *   edgedata: cell array of edge data
 *   schedule: array of task structs where each
 *     task struct:
 *       vertex: numeric vertex id
 *       update function: string. Name of the update function
 * 
 *  Returns 0 on     
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  updates_initialize();
  
  
  
  
}
