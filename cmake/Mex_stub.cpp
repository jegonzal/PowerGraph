// Adopted from: http://www.cmake.org/Wiki/images/7/72/Mex_stub.cpp
// and tutorial: http://www.cmake.org/Wiki/CMake:MatlabMex
// on June 9, 2010  (akyrola)

#include "mex.h"


extern void __mexFunction__(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
extern void __at_exit__();

static void at_exit();

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    
  mexAtExit(&at_exit);

  __mexFunction__(nlhs, plhs, nrhs, prhs);

}

static void at_exit()
{
  __at_exit__();    
}
