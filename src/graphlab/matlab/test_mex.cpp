#include <mex.h>
#include "rtwtypes.h"
#include "updates_types.h"
#include "mx_emx_converters.hpp"
#include "updates_initialize.h"
#include "gl_emx_graphtypes.hpp"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  gl_emx_vertextype v;
  gl_emx_edgetype e;
  if (nlhs != 1 || nrhs != 1) {
    mexPrintf("args invalid\n");
    return;
  }
  updates_initialize();
}
