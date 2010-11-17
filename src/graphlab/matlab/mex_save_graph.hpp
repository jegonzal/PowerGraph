#ifndef GRAPHLAB_MEX_HPP
#define GRAPHLAB_MEX_HPP

struct mex_parameters{
  const mxArray* vdata;
  const mxArray* adjmat;
  const mxArray* edata;
  const mxArray* options;
  const mxArray* graphfile;
  const mxArray* strict;
};



#endif