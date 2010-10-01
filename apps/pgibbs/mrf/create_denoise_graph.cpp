
/**
 *
 * Create a noisy instance to solve
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <set> 
#include <algorithm>
#include <limits>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/bind.hpp>

#include <graphlab.hpp>

// Image reading/writing code
#include "image.hpp"
#include "data_structures.hpp"
#include "image_denoise.hpp"



// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>


// Command Line Parsing =======================================================>
struct options {
  size_t num_rings;
  size_t rows;
  size_t cols;
  double sigma;
  double lambda;  
  std::string smoothing;
  std::string drawing;
  std::string corruption;
};


/**
 * Parse the command line arguments.  Returns false if there was a
 * problem in parsing command line arguments
 */    
bool parse_command_line(int argc, char** argv, options& opts);

/**
 * Display the program options
 */
void display_options(options& opts);

// MAIN =======================================================================>
int main(int argc, char** argv) {
  
  std::cout << "This program uses gibbs sampling to denoise a synthetic image."
            << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

  // Parse command line arguments --------------------------------------------->
  options opts; 
  bool success = parse_command_line(argc, argv, opts);
  if(!success)  return EXIT_FAILURE;
  display_options(opts);


  // Create a denoise problem struct
  denoise_problem problem(opts.num_rings,
                          opts.rows,
                          opts.cols,
                          opts.sigma,
                          opts.lambda,
                          opts.smoothing,
                          opts.drawing,
                          opts.corruption);

  problem.save_all();

  std::cout << problem << std::endl;

    
  
  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;


  
} // End of main





bool parse_command_line(int argc, char** argv, options& opts) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Denoise a randomly generated image using Gibbs Sampling.");
  // Set the program options
  desc.add_options()
    ("help",   "produce this help message")
    ("rings",  boost_po::value<size_t>(&(opts.num_rings))->default_value(5),
     "Number of rings in the noisy image")
    ("rows",  boost_po::value<size_t>(&(opts.rows))->default_value(200),
     "Number of rows in the noisy image")
    ("cols",  boost_po::value<size_t>(&(opts.cols))->default_value(200),
     "Number of columns in the noisy image")
    ("sigma",  boost_po::value<double>(&(opts.sigma))->default_value(1.2),
     "Standard deviation of noise.")
    ("lambda",  boost_po::value<double>(&(opts.lambda))->default_value(3),
     "Smoothness parameter (larger => smoother).")    
    ("smoothing",
     boost_po::value<std::string>(&(opts.smoothing))->default_value("square"),
     "Options are {square, laplace}")
    ("corruption",
     boost_po::value<std::string>(&(opts.corruption))->default_value("gaussian"),
     "{gaussian, flip}")
    ("drawing",
     boost_po::value<std::string>(&(opts.drawing))->default_value("sunset"),
     "Options are {sunset, checkerboard, ising}");

  // Parse the arguments
  boost_po::variables_map vm;
  boost_po::positional_options_description pos_opts;
  store(boost_po::command_line_parser(argc, argv)
        .options(desc).positional(pos_opts).run(), vm);
  boost_po::notify(vm);
  if(vm.count("help") ) {
    std::cout << desc << std::endl;
    return false;
  }
  return true;
} // end of parse command line arguments


void display_options(options& opts) {
  std::cout << "num_rings:      " << opts.num_rings << std::endl
            << "rows:           " << opts.rows << std::endl
            << "cols:           " << opts.cols << std::endl
            << "sigma:          " << opts.sigma << std::endl
            << "lambda:         " << opts.lambda << std::endl
            << "smoothing:      " << opts.smoothing << std::endl
            << "drawing:        " << opts.drawing << std::endl
            << "corruption:     " << opts.corruption << std::endl;
}

