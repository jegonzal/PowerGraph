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
 */


#ifndef _CG_HPP
#define _CG_HPP

#include <graphlab.hpp>
#include "math.hpp"
#include <graphlab/macros_def.hpp>

extern uint32_t m,n;
extern gl_types::core * glcore;
extern bool debug;
extern bool square;
extern int cg_maxiter;
extern bool cg_resid;
extern bool cg_noop;

void init_row_cols();

using namespace graphlab;

/*function [x] = conjgrad(A,b,x)
    r=b-A*x; ///DIST
    p=r;     //SER
    rsold=r'*r;  //SER

    for i=1:size(A,1)
        Ap=A*p;               //DIST
        alpha=rsold/(p'*Ap);  //SER
        x=x+alpha*p;          //SER
        r=r-alpha*Ap;         //SER
        rsnew=r'*r;           //SER
        if sqrt(rsnew)<1e-10  //SER
              break;
        end
        p=r+rsnew/rsold*p;    //SER
        rsold=rsnew;          //SER
    end
end
*/

double cg(gl_types::core * _glcore, std::vector<double> & means, double &diff){

    glcore = _glcore;
    int tmp=n; n=m; m = tmp;
    diff = NAN;
    init_row_cols();

    DistMat A;
    DistVec b(0,true), prec(1,true), r(2, true), p(3,true), x(4,true), Ap(5), t(6,true);
    //initialize startng guess
    if (!square)
      x = ones(0.5,m);
    else
      x = ones(0.5,n);
    
    DistDouble rsold, rnew, alpha;

    /* r = -A*x+b;
       if (~sqaure)
         r=A'*r;
       end
       p = r;
       rsold = r'*r;
    */
    r=-A*x+b; 
    if (!square)
      r=A._transpose()*r;
    p = r;
    rsold = r._transpose()*r;

     /*
     for i=1:size(A,1)
        Ap=A*p; 
        if (~square)
          Ap=A'*Ap; 
        end             
        alpha=rsold/(p'*Ap); 
        x=x+alpha*p;        
        r=r-alpha*Ap;      
        rsnew=r'*r;       
        if sqrt(rsnew)<1e-10 
              break;
        end
        p=r+rsnew/rsold*p;  
        rsold=rsnew;       
    end
    */

    for (int i=1; i <= std::min(cg_maxiter,size(A,1)); i++){
        Ap=A*p;
        if (!square)
          Ap= A._transpose()*Ap;
        alpha = rsold/(p._transpose()*Ap);
        x=x+alpha*p;
   
        if (cg_resid){
          t=A*x-b;
          logstream(LOG_INFO)<<"Iteration " << i << " approximated solution redidual is " << norm(t).toDouble() << std::endl;
        }

        r=r-alpha*Ap;
        rnew=r._transpose()*r;
        if (sqrt(rnew)<1e-10){
          diff = sqrt(rnew).toDouble();
          logstream(LOG_INFO)<<" Conjugate gradient converged in iteration "<<i<<" to an accuracy of "  << diff << std::endl; 
          break;
        }
        p=r+(rnew/rsold)*p;
        rsold=rnew;
    }

    x.to_vec(means);
    return runtime;

}
#include <graphlab/macros_undef.hpp>
#endif
