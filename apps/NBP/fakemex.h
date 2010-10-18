#ifndef FAKE_MEX
#define FAKE_MEX

#include <stdlib.h>

#define mxMalloc malloc
#define mxFree free
#define mexErrMsgTxt printf
#define mxDestroyArray

#endif
