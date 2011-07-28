/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
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
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *  
 *  Code written by Danny Bickson, CMU
 *  Any changes to the code must include this original license notice in full.
 *  This file implements the shooting algorithm for solving Lasso problem
 */


#ifndef _LASSO_HPP
#define _LASSO_HPP

extern bool debug;

vec lasso(mat A, vec b, double lambda, int max_iter, int D){

  assert(lambda > 0);
  vec beta;

  assert(A.rows() == A.cols());
  assert(A.rows() == D);
  assert(b.size() == D);
  bool ret = itpp::ls_solve(A, b, beta);
  assert(ret);

  if (debug)
     cout<<"initial beta is: " << beta << endl;

  for (int i=0; i< max_iter; i++){
	for (int j=0; j<D; j++){
            vec proj = A*beta;
            if (debug)
                cout<< " proj[j] " << proj[j] << " 2*beta[j] " << 2*beta[j] << " b[j] " << b[j] << endl;
            double S0 = proj[j] - 2*beta[j] - b[j];
            if (S0 > lambda)
                beta[j] = (lambda-S0)/2;
	    else if (S0 < -lambda)
		beta[j] = (-lambda -S0)/2;
            else
		beta[j] = 0;

	    if (debug)
               cout<<" round " << i << " dim " << j << " S0: " << S0 << " beta: " << beta << endl;

        }
  }
  if (sum(abs(beta)) == 0)
    logstream(LOG_ERROR) << " --lasso_lambda=XX value is too large, resulting in a zero solution. Try to decrease value and try again" << endl;
  return beta;

}







#endif
