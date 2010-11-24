/**
 * This file contains an example of graphlab used for discrete loopy
 * belief propagation in a pairwise markov random field to denoise a
 * synthetic noisy image.
 *
 *  \author Yucheng
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <iostream>
#include <map>

#include <graphlab.hpp>
#include <limits>
#include <float.h>
#include "image.hpp"

// include itpp
#include <itpp/itstat.h>
#include <itpp/itbase.h>

#include "gaussian_mixture.hpp"
// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

#define PROPOSAL_STDEV 20

using namespace itpp;
using namespace std;



// MAIN =======================================================================>
int main(int argc, char** argv) {
  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  size_t iterations = 100;
  size_t numparticles = 100;
  std::string pred_type = "map";

  std::string gmmfile= "";
  std::string inputfile = "";


  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Loopy BP image denoising");
  clopts.attach_option("iterations",
                       &iterations, iterations,
                       "Number of iterations");
  clopts.attach_option("particles",
                       &numparticles, numparticles,
                       "Number of particlesw");
  clopts.attach_option("gmmfile",
                       &gmmfile, std::string(""),
                       "gmm mixture file");
  clopts.attach_option("inputfile",
                       &inputfile, std::string(""),
                       "input file");
  clopts.attach_option("damping",
                       &damping, damping,
                       "damping");
  clopts.attach_option("resample",
                       &RESAMPLE_FREQUENCY, RESAMPLE_FREQUENCY,
                       "resampling frequency");
  clopts.attach_option("mcmc",
                       &MCMCSTEPS, MCMCSTEPS,
                       "mcmc steps per resample");

  // set default scheduler type
  clopts.scheduler_type = "splash(100)";
  clopts.scope_type = "edge";

  bool success = clopts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }

  // fill the global vars
  MAX_ITERATIONS = iterations;

  // load the potentials mixture components
  it_ifile f(gmmfile.c_str());

  // weigghts
  mat edgecenter, edgesigma, edgeweight;
  mat nodecenter, nodesigma, nodeweight;
  ivec truedata;
  ivec imgsize;
  // intermediate types to use...
  imat integermat;
  vec doublevec;
  f >> Name("edge_ce") >> integermat;   edgecenter = to_mat(integermat);
  f >> Name("edge_alpha") >> doublevec; edgeweight = doublevec;
  f >> Name("edge_sigma") >> doublevec; edgesigma = doublevec;

  f >> Name("like_ce") >> nodecenter;
  f >> Name("like_alpha") >> doublevec; nodeweight = doublevec;
  f >> Name("like_sigma") >> doublevec; nodesigma = doublevec;
  f >> Name("img1") >> truedata;
  f >> Name("isize") >> imgsize;

  size_t rows = imgsize(0);
  size_t cols = imgsize(1);
  std::cout << "Image size is "
            << rows << " x " << cols << std::endl;
  // make the GMMs
  GaussianMixture<2> edgepot(edgecenter, edgesigma, edgeweight);
  GaussianMixture<2> nodepot(nodecenter, nodesigma, nodeweight);
  edgepot.simplify();
  edgepot.print();
  // convert the true image to an image
  image trueimg(rows, cols);
  for (size_t i = 0;i < truedata.size(); ++i) {
    trueimg.pixel(i) = truedata(i);
  }

  // load the observations
  it_ifile imgfile(inputfile.c_str());
  vec observations;
  imgfile >> Name("obs2") >> observations;
  // convert observations to an image
  image img(rows, cols);
  for (size_t i = 0;i < observations.size(); ++i) {
    img.pixel(i) = observations(i);
  }
  img.save("noisy.pgm");
  trueimg.save("source_img.pgm");

  // Create the graph --------------------------------------------------------->
  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);

  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;

  construct_graph(img, edgepot, nodepot, numparticles, core.graph());



  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function", bp_update);

  std::cout << "Running the engine. " << std::endl;


  // Add the bp update to all vertices
  core.add_task_to_all(bp_update, 100.0);
  // Starte the engine
  double runtime = core.start();

  size_t update_count = core.last_update_count();
  std::cout << "Finished Running engine in " << runtime
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;


  // Saving the output -------------------------------------------------------->
  std::cout << "Rendering the cleaned image. " << std::endl;

  for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    float a = vdata.max_asg();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = size_t(a);
  }
  img.save("pred_map.pgm");


  for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    float a = vdata.average();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = (a);
  }
  img.save("pred_exp.pgm");


  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
} // End of main

