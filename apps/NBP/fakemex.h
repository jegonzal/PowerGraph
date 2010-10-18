#ifndef FAKE_MEX
#define FAKE_MEX

#include <stdlib.h>
#include <assert.h>

#define mxMalloc malloc
#define mxFree free
#define mexErrMsgTxt(a) {printf(a); assert(false); }
#define mxDestroyArray

#endif
