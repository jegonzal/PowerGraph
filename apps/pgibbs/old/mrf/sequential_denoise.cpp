/**
 *
 * This file contains an experimental version of blocked gibbs
 * denoising.  The current implementation does not use the GraphLab
 * engine but relies on the graphlab support code.
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
#include "drawing.hpp"
#include "sequential_gibbs.hpp"
#include "sequential_tree_gibbs.hpp"
#include "parallel_tree_gibbs.hpp"


// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>




void compare_samplers(denoise_problem& problem,
                      size_t max_tree_size,
                      size_t step_count,
                      size_t step_size) {
  // Make two copies of the graph
  std::cout << "Copying graphs" << std::endl;
  gl::graph sweep_graph = problem.graph;
  gl::graph tree_graph = problem.graph;

  graphlab::unary_factor marginal;
  size_t nsamples = 0;
  for(size_t step = 0; step < step_count; ++step) {
    nsamples += step_size;
    std::cout << "------------------------------------------------------"
              << std::endl;
    std::cout << "Nsamples = " << nsamples << std::endl;
 
    std::cout << "Running forward sweep" << std::endl;
    random_sweep(sweep_graph, problem.edge_factor, step_size);
    //forward_sweep(sweep_graph, problem.edge_factor, step_size);
    
    std::cout << "Running tree sampler" << std::endl;
    tree_sample(tree_graph, problem.edge_factor, max_tree_size, step_size);
    // comb_sample(tree_graph, rows, cols, edge_factor);

    std::cout << "Saving Images" << std::endl;
    image sweep_img(problem.rows, problem.cols);
    image tree_img(problem.rows, problem.cols);
    image pair_img(problem.rows, problem.cols*2);
    image tree_updates_img(problem.rows, problem.cols);
    for(size_t v = 0; v < problem.graph.num_vertices(); ++v) {
      std::pair<size_t, size_t> loc = sweep_img.loc(v);
      double pixel_value = 0;

      marginal = sweep_graph.vertex_data(v).belief;
      marginal.normalize();
      pixel_value = marginal.expectation();
      //  pixel_value = sweep_graph.vertex_data(v).asg;
      sweep_img.pixel(v) = pixel_value;
      pair_img.pixel(loc.first, loc.second) = pixel_value;

      marginal = tree_graph.vertex_data(v).belief;
      marginal.normalize();
      pixel_value = marginal.expectation();
      //  pixel_value = tree_graph.vertex_data(v).asg;
      tree_img.pixel(v) = pixel_value;
      pair_img.pixel(loc.first, loc.second + problem.cols) = pixel_value;
      tree_updates_img.pixel(v) = tree_graph.vertex_data(v).updates;
    }
    sweep_img.save(make_filename("sweep", ".pgm", nsamples).c_str());
    tree_img.save(make_filename("tree", ".pgm", nsamples).c_str());
    pair_img.save(make_filename("pair", ".pgm", nsamples).c_str());
    tree_updates_img.save(make_filename("tree_updates", ".pgm", nsamples).c_str());
    std::cout << tree_updates_img.min() << ", "
              << tree_updates_img.max() << std::endl;
  }


  save_beliefs(sweep_graph, "sweep_blfs.tsv");
  save_beliefs(tree_graph,  "tree_blfs.tsv");

  
} // end of compare samplers






// Command Line Parsing =======================================================>
struct options {
  size_t samples;
  size_t num_rings;
  size_t rows;
  size_t cols;
  double sigma;
  double lambda;

  size_t max_tree_size;
  size_t step_count;

  std::string problem;
  
  std::string smoothing;
  std::string orig_fn;
  std::string noisy_fn;
  std::string pred_fn;
  std::string pred_type;
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

  // Seed The random number generator 
  graphlab::random::seed();
  
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

  denoise_problem problem;
  if(opts.problem == "") {
    problem = denoise_problem(opts.num_rings,
                              opts.rows,
                              opts.cols,
                              opts.sigma,
                              opts.lambda,
                              opts.smoothing);
  } else {
    problem.load(opts.problem);
  }

  std::cout << problem << std::endl;
  
  
  // Create synthetic images -------------------------------------------------->
  // Creating image for denoising
  std::cout << "Creating a synethic image." << std::endl;

  std::cout << "Saving image. " << std::endl;
  problem.original.save(opts.orig_fn.c_str());

  std::cout << "Saving corrupted image. " << std::endl;
  problem.noisy.save(opts.noisy_fn.c_str());

  

  // std::cout << "Running the gibbs sampler forever" << std::endl;
  // gl::graph true_graph = graph;
  // forward_sweep(true_graph, edge_potential, true_graph.num_vertices() * 10000);
  // draw_graph(true_graph, "truth.pgm", opts.rows, opts.cols);
  // save_beliefs(true_graph, "true_blfs.tsv");
  
  
  //   std::cout << "Drawing some trees" << std::endl;
  // draw_trees(graph, opts.max_tree_size, 20, img.rows(), img.cols());



  // std::vector<vertex_id_t> tree_order;
  // grow_comb_tree(graph, false, opts.rows, opts.cols, tree_order);
  // draw_tree(tree_order, opts.rows, opts.cols, "combtree_1.pgm");
  // grow_comb_tree(graph, true, opts.rows, opts.cols, tree_order);
  // draw_tree(tree_order, opts.rows, opts.cols, "combtree_2.pgm");


  size_t num_verts = problem.graph.num_vertices();

  compare_samplers(problem,
                   opts.max_tree_size,
                   opts.step_count,
                   5 * num_verts);




  // slow_compare_samplers(graph,
  //                  edge_factor,
  //                  opts.rows, opts.cols,
  //                  num_verts * 50,
  //                  num_verts * 5);



  
  // Run the sampler ----------------------------------------------------------->

//   for(size_t i = 0; i < opts.samples; ++i) {
//     forward_sweep(graph, edge_factor);
//   }
  // tree_sample(graph, edge_potential, opts.samples);
  // for(size_t i = 0; i < graph.num_vertices(); ++i) 
  //   graph.vertex_data(i).belief.normalize();



  

  // // Create the edge weights -------------------------------------------------->
  // std::cout << "Computing edge weights. " << std::endl;
  // update_edge_weights(graph, edge_potential);
  // std::cout << "Saving edge weights. " << std::endl;
  // save_edge_weights(graph, img.rows(), img.cols());
  // std::cout << "Finished!" << std::endl;

  

  // // Saving the output -------------------------------------------------------->
  // std::cout << "Rendering the cleaned image. " << std::endl;
  // if(opts.pred_type == "map") {
  //   for(size_t v = 0; v < graph.num_vertices(); ++v) {
  //     const vertex_data& vdata = graph.vertex_data(v);
  //     img.pixel(v) = vdata.belief.max_asg();    
  //   }
  // } else if(opts.pred_type == "exp") {
  //   for(size_t v = 0; v < graph.num_vertices(); ++v) {
  //     const vertex_data& vdata = graph.vertex_data(v);    
  //     img.pixel(v) = vdata.belief.expectation();
  //   }
  // } else {
  //   std::cout << "Invalid prediction type! : " << opts.pred_type
  //             << std::endl;
  //   return EXIT_FAILURE;
  // }
  // std::cout << "Saving cleaned image. " << std::endl;
  // img.save(opts.pred_fn.c_str());
  
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
    ("samples",  boost_po::value<size_t>(&(opts.samples))->default_value(100),
     "Number of samples to generate")
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
    
    ("treesize", boost_po::value<size_t>(&(opts.max_tree_size))->
     default_value(200*200),
     "Maximum tree size for tree sampling")

    ("stepcount", boost_po::value<size_t>(&(opts.step_count))->
     default_value(50),
     "number of sampling steps")
    
    ("problem",
     boost_po::value<std::string>(&(opts.problem))->default_value(""),
     "A problem description file")


    
    ("smoothing",
     boost_po::value<std::string>(&(opts.smoothing))->default_value("square"),
     "Options are {square, laplace}")

    ("orig",
     boost_po::value<std::string>(&(opts.orig_fn))->default_value("source_img.pgm"),
     "Original image file name.")
    ("noisy",
     boost_po::value<std::string>(&(opts.noisy_fn))->default_value("noisy_img.pgm"),
     "Noisy image file name.")
    ("pred",
     boost_po::value<std::string>(&(opts.pred_fn))->default_value("pred_img.pgm"),
     "Predicted image file name.")
    ("pred_type",
     boost_po::value<std::string>(&(opts.pred_type))->default_value("map"),
     "Predicted image type {map, exp}");
  // Parse the arguments
  boost_po::variables_map vm;
  boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
  boost_po::notify(vm);
  if(vm.count("help")) {
    std::cout << desc << std::endl;
    return false;
  }
  return true;
} // end of parse command line arguments


void display_options(options& opts) {
  std::cout << "samples:        " << opts.samples << std::endl
            << "num_rings:      " << opts.num_rings << std::endl
            << "rows:           " << opts.rows << std::endl
            << "cols:           " << opts.cols << std::endl
            << "sigma:          " << opts.sigma << std::endl
            << "lambda:         " << opts.lambda << std::endl
            << "treesize:       " << opts.max_tree_size << std::endl
            << "step_count:     " << opts.step_count << std::endl
            << "problem:        " << opts.problem << std::endl
            << "smoothing:      " << opts.smoothing << std::endl           
            << "orig_fn:        " << opts.orig_fn << std::endl
            << "noisy_fn:       " << opts.noisy_fn << std::endl
            << "pred_fn:        " << opts.pred_fn << std::endl
            << "pred_type:      " << opts.pred_type << std::endl;
}




