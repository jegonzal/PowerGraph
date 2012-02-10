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
bool debug = false;
int max_iter = 10;

struct vertex_data {

  vec pvec;
  double A_ii;
  vertex_data() : A_ii(1) {
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

int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile, yfile, xfile;
  std::string format = "matrixmarket";
  real_type threshold = 1e-5;
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

  //unit testing
  if (unittest == 1){
/*A= [  1.8147    0.9134    0.2785
       0.9058    1.6324    0.5469
       0.1270    0.0975    1.9575 ];

  y= [ 0.9649    0.1576    0.9706 ]';
  x= [ 0.6803   -0.4396    0.4736 ]';
*/
    datafile = "A"; yfile = "y"; xfile = "x"; sync_interval = 120;
  }

  std::cout << "Load matrix A" << std::endl;
  bipartite_graph_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());
  std::cout << "Load Y values" << std::endl;
  load_vector(yfile, format, matrix_info, core.graph(), CG_X, false);
  std::cout << "Load x values" << std::endl;
  load_vector(xfile, format, matrix_info, core.graph(), CG_Y, true);
  
  math_info mi;
  init_math(&core.graph(), &core);
  //aggregator acum;
  //core.add_aggregator("sync", acum, sync_interval);
  //core.add_global("REAL_NORM", double(0));
  //core.add_global("RELATIVE_NORM", double(0));
  //core.add_global("THRESHOLD", threshold); 

  
  //  int tmp=ps.n; ps.n=ps.m; ps.m = tmp;
  //  diff = NAN;
  //  init_row_cols();

    DistMat A(matrix_info,mi);
    DistVec b(mi, matrix_info, CG_Y,true, "b");
    DistVec prec(mi, matrix_info, CG_PREC,true, "prec");
    DistVec r(mi, matrix_info, CG_R, true, "r");
    DistVec p(mi, matrix_info, CG_P,true, "p");
    DistVec x(mi, matrix_info, CG_X,true, "x");
    DistVec Ap(mi, matrix_info, CG_AP, false, "Ap");
    DistVec t(mi, matrix_info, CG_T,true, "t");
    //initialize startng guess

    vec init_x;
    if (debug)
      init_x = ones(matrix_info.num_nodes(false));
    x = init_x;

    DistDouble rsold(mi), rnew(mi), alpha(mi);

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

  
  vec ret = fill_output(&core.graph(), matrix_info, CG_X);
  write_output_vector(datafile + "x.out", format, ret, false);


  if (unittest == 1){
    double real_norm = core.get_global<double>("REAL_NORM");
    double relative_norm = core.get_global<double>("RELATIVE_NORM");
    assert(real_norm < 1e-30);
    assert(relative_norm < 1e-30);
  }

   return EXIT_SUCCESS;
}



