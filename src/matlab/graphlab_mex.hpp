#ifndef GRAPHLAB_MEX_HPP
#define GRAPHLAB_MEX_HPP
#include "mex.h"
#include <graphlab.hpp>

typedef graphlab::graph<mxArray*, mxArray*> matlab_graph;
typedef graphlab::types<matlab_graph> gl;

typedef int32_t (*matlab_update_function)(int32_t, int32_t);

#endif