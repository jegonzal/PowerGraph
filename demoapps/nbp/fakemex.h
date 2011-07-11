#ifndef FAKE_MEX
#define FAKE_MEX

/* Written by Danny Bickson, CMU */
/* File for defining mex commands in their C implementation */


#include <stdlib.h>
#include <assert.h>

#define mxMalloc malloc
#define mxFree free
#define mexErrMsgTxt(a) {printf(a); assert(false); }
#define mxDestroyArray

#endif
