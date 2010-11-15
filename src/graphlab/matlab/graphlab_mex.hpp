#ifndef GRAPHLAB_MEX_HPP
#define GRAPHLAB_MEX_HPP
#include <mex.h>

struct mex_parameters{
  const mxArray* vdata;
  const mxArray* adjmat;
  const mxArray* edata;
  const mxArray* sched;
  const mxArray* strict;
};



#endif