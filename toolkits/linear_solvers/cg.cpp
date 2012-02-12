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



#include <graphlab.hpp>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include "../shared/mathlayer.hpp"

using namespace graphlab;
double regularization = 0;
bool calc_residual = false;
bool final_residual = true;
double threshold = 1e-5;
bool debug = false;
int max_iter = 10;
int data_size = 7;
std::string datafile, yfile, xfile;

struct vertex_data {

  vec pvec;
  double A_ii;
  vertex_data() : A_ii(1) {
    pvec = zeros(data_size);
  }
  void add_self_edge(double value) { A_ii = value + regularization; }

  void set_val(double value, int field_type) { 
    pvec[field_type] = value;
  }
  //only one output for jacobi - solution x
  double get_output(int field_type){ return pvec[field_type]; }
}; // end of vertex_data

struct edge_data {
  real_type weight;
  edge_data(double weight = 0) : weight(weight) { }
};


typedef graphlab::graph<vertex_data, edge_data> graph_type;
#include "../shared/math.hpp"

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
enum input_pos{
  CG_Y = 0,
  CG_PREC = 1,
  CG_R = 2,
  CG_P = 3,
  CG_X = 4,
  CG_AP = 5,
  CG_T = 6
};


void test_math(int unittest, bipartite_graph_descriptor & info, math_info & mi){
    DistMat A(info);
    DistVec b(info, CG_Y,true, "b");
    DistVec prec(info, CG_PREC,true, "prec");
    DistVec p(info, CG_P,true, "p");
    DistVec x(info, CG_X,true, "x");
    DistVec Ap(info, CG_AP, false, "Ap");
 
    x = ones(3);
    assert(x[0] == 1);
    assert(x[1] == 1);
    assert(x[2] == 1);
    x = -x;
    assert(x[0] == -1);
    assert(x[1] == -1);
    assert(x[2] == -1);
    DistDouble factor(2);
    assert(factor == 2);
    x = factor*x;
    assert(x[0] == -2);
    assert(x[1] == -2);
    assert(x[2] == -2);
    x = x + b;
    vec bret = b.to_vec();
    vec ret = x.to_vec();
    assert(ret[0] == bret[0] - 2);
    assert(ret[1] == bret[1] - 2);
    assert(ret[2] == bret[2] - 2);

    p = A*b;
    vec pret = p.to_vec();
    assert(pow(p[0] - 0.89646579890448475,2) <1e-15); 

    Ap = A._transpose() * p;
    assert(pow(Ap[0] -1.69277,2)<1e-8); 

    x = zeros(3);
    p = A*x;
    assert(x[0] == 0);
    assert(x[2] == 0);
    assert(p[2] == 0);
    assert(p[0] == 0);

    x = zeros(3);
    assert(x.size() == 3);
    p = A*x + b;
    assert(p.size() == 3);
    //0.680272 -0.439596 0.473605
    assert(pow(p[0] - 0.680272,2)<1e-8);
    assert(pow(p[1] - -0.439596,2)<1e-8);
    assert(pow(p[2] - 0.473605,2)<1e-8);

    p = b + A*x;
    assert(p.size() == 3);
    assert(pow(p[0] - 0.680272,2)<1e-8);
    assert(pow(p[1] - -0.439596,2)<1e-8);
    assert(pow(p[2] - 0.473605,2)<1e-8);

    p = factor*b+factor*A*x;
    assert(p.size() == 3);
    assert(pow(p[0] - (2*0.680272),2)<1e-8);
    assert(pow(p[1] - (2*-0.439596),2)<1e-8);
    assert(pow(p[2] - (2*0.473605),2)<1e-8);

    p = b*factor + A*x*factor;
    assert(p.size() == 3);
    assert(pow(p[0] - (2*0.680272),2)<1e-8);
    assert(pow(p[1] - (2*-0.439596),2)<1e-8);
    assert(pow(p[2] - (2*0.473605),2)<1e-8);

    exit(0);
}

void setup_unittest(int unittest){   
   if (unittest == 1){
     datafile = "unittest1.mtx"; yfile = "unittest1.mtxv"; max_iter = 20; threshold = 1e-10;
   }
}

void verify_values(int unittest, double residual){
   if (unittest == 1)
     assert(residual < 1e-14);
}

int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string format = "matrixmarket";
  size_t sync_interval = 10000;
  int unittest = 0;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("yfile", &yfile, yfile,
                       "vector y input file");
  clopts.attach_option("xfile", &xfile, xfile,
                       "vector x input file (optional)");
  clopts.attach_option("threshold", &threshold, threshold, "termination threshold.");
  clopts.add_positional("threshold");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("syncinterval", 
                       &sync_interval, sync_interval, 
                       "sync interval (number of update functions before convergen detection");
  clopts.attach_option("regularization", &regularization, regularization, 
	               "regularization added to the main diagonal");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("max_iter", &max_iter, max_iter, "max number of iterations");
  clopts.attach_option("calc_residual", &calc_residual, calc_residual, "calc residual in each iteration"); 
  clopts.attach_option("final_residual", &final_residual, final_residual, "calc residual at the end (norm(Ax-b))");
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab Linear solver library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Gaussian Belief Propagation, "
    << "Jacobi method, Conjugate Gradient" << std::endl;



  // Create a core
  graphlab::core<graph_type, Axb> core;
  core.set_options(clopts); // Set the engine options

  if (unittest > 0)
    setup_unittest(unittest);

  std::cout << "Load matrix A" << std::endl;
  bipartite_graph_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());
  std::cout << "Load Y values" << std::endl;
  load_vector(yfile, format, matrix_info, core.graph(), CG_X, false);
  std::cout << "Load x values" << std::endl;
  load_vector(xfile, format, matrix_info, core.graph(), CG_Y, true);
  
  math_info mi;
  init_math(&core.graph(), &core);

  if (unittest == 2)
    test_math(unittest, matrix_info, mi); 

    DistMat A(matrix_info);
    DistVec b(matrix_info, CG_Y,true, "b");
    DistVec prec(matrix_info, CG_PREC,true, "prec");
    DistVec r(matrix_info, CG_R, true, "r");
    DistVec p(matrix_info, CG_P,true, "p");
    DistVec x(matrix_info, CG_X,true, "x");
    DistVec Ap(matrix_info, CG_AP, false, "Ap");
    DistVec t(matrix_info, CG_T,true, "t");
    //initialize startng guess

    vec init_x;
    if (debug)
      init_x = ones(matrix_info.num_nodes(false));
    else
      init_x = randu(matrix_info.num_nodes(false));

    x = init_x;

    DistDouble rsold, rnew, alpha;

    /* r = -A*x+b;
       if (~sqaure)
         r=A'*r;
       end
       p = r;
       rsold = r'*r;
    */
    r=-A*x+b; 
    if (!matrix_info.is_square())
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

    for (int i=1; i <= std::min(max_iter,size(A,1)); i++){
        Ap=A*p;
        if (!matrix_info.is_square())
          Ap= A._transpose()*Ap;
        alpha = rsold/(p._transpose()*Ap);
        x=x+alpha*p;
   
        if (calc_residual){
          t=A*x-b;
          logstream(LOG_INFO)<<"Iteration " << i << " approximated solution residual is " << norm(t).toDouble() << std::endl;
        }

        r=r-alpha*Ap;
        rnew=r._transpose()*r;
        if (sqrt(rnew)< threshold){
          double diff = sqrt(rnew).toDouble();
          logstream(LOG_INFO)<<" Conjugate gradient converged in iteration "<<i<<" to an accuracy of "  << diff << std::endl; 
          break;
        }
        p=r+(rnew/rsold)*p;
        rsold=rnew;
    }

  if (final_residual){
    p = A*x -b;
    DistDouble ret = norm(p);
    logstream(LOG_INFO) << "Solution converged to residual: " << ret.toDouble() << std::endl;
    if (unittest > 0){
     verify_values(unittest, ret.toDouble());
    }
  }
 

  vec ret = fill_output(&core.graph(), matrix_info, CG_X);
  write_output_vector(datafile + "x.out", format, ret, false);


   return EXIT_SUCCESS;
}



