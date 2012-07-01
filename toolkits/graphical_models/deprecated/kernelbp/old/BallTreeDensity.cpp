/* Copyright (c) 2003 Alexander Ihler
 * Original code from: http://www.ics.uci.edu/~ihler/code/index.html
 *
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License. */

 //
// Matlab MEX interface for KD-tree C++ functions
//
// Written by Alex Ihler and Mike Mandel
// Copyright (C) 2003 Alexander Ihler
//

//#define MEX
#include "cpp/BallTreeDensity.h"
#ifdef MEX
#include "mex.h"
#endif

#ifdef MEX
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // check for the right number of arguments
  if((nrhs < 3)||(nrhs > 4))
    mexErrMsgTxt("Takes 3-4 input arguments");
  if(nlhs != 1)
    mexErrMsgTxt("Outputs one result (a structure)");

  if (nrhs == 3) //                          points, weights, bandwidths
    plhs[0] = BallTreeDensity::createInMatlab(prhs[0],prhs[1],prhs[2]);
  else {          //                          points, weights, bandwidths,type
    int ktype = (int) mxGetScalar(prhs[3]);
    plhs[0] = BallTreeDensity::createInMatlab(prhs[0],prhs[1],prhs[2],(BallTreeDensity::KernelType) ktype);
  }
}
#else


#endif
