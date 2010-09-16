
/**
 * This runs the parallel blocked gibbs samplering using graphlab
 */

#include <iostream>


#include <graphlab.hpp>



#include "data_structures.hpp"




int main(int argc, char** argv) {
  std::cout << "Running parallel blocked gibbs on a uniform ising model."
            << std::endl;


  size_t rows = 20;
  size_t cols = 20;
  size_t colors = 2;
  double lambda = 1;

  size_t nsamples = 100;



  // Setup the command line arguments processor
  graphlab::command_line_options clopts("PGibbs Ising Solver");
  clopts.attach_option("rows", 
                       &rows, rows,
                       "Number of rows in grid");
  clopts.attach_option("cols", 
                       &cols, cols,
                       "Number of cols in grid");
  clopts.attach_option("colors", 
                       &colors, colors,
                       "Number of colors");

  clopts.attach_option("lambda", 
                       &lambda, lambda,
                       "Edge parameter f(a,b) = exp(-lambda * I[a == b])");
  
  clopts.attach_option("nsamples",
                       &nsamples, nsamples,
                       "Total number of samples = nsamples * rows * cols");

  
  // Set the default scheduler to be fifo
  clopts.scheduler_type = "fifo";
  // Set the default scope to be edge
  clopts.scope_type = "edge";

  // Parse the command line arguments
  if( !clopts.parse(argc, argv) ) {
    std::cout << "Error in parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }

  
  // Create a graphlab core initialized with the command line settings
  gl::core core;
  core.set_engine_options(clopts);





  return EXIT_SUCCESS;
}
