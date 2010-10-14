/**
 * Run parallel junction tree gibbs sampling on a factorized model
 */

#include <cstdlib>
#include <iostream>


#include <graphlab.hpp>


#include "image.hpp"
#include "data_structures.hpp"
#include "sequential_jt_gibbs.hpp"



#include "jt_worker.hpp"





#include <graphlab/macros_def.hpp>




int main(int argc, char** argv) {
  std::cout << "This program runs junction tree blocked MCMC "
            << "inference on large factorized models."
            << std::endl;

  std::string model_filename = "";

  size_t treesize = 1000;
  bool priorities = false;
  float runtime = 10;

  // Command line parsing
  graphlab::command_line_options clopts("Parallel Junction Tree MCMC");
  clopts.attach_option("model", 
                       &model_filename, model_filename,
                       "Alchemy formatted model file");
  clopts.add_positional("model");

  clopts.attach_option("runtime", 
                       &runtime, runtime,
                       "total runtime in seconds");

  clopts.attach_option("treesize", 
                       &treesize, treesize,
                       "The number of variables in a junction tree");
  clopts.attach_option("priorities",
                       &priorities, priorities,
                       "Use priorities?");



  clopts.scheduler_type = "fifo";
  clopts.scope_type = "edge";
  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Load alchemy file." << std::endl;
  factorized_model factor_graph;
  factor_graph.load_alchemy(model_filename);

  std::cout << "Building graphlab MRF." << std::endl;
  mrf::graph_type mrf_graph;
  construct_mrf(factor_graph, mrf_graph);
  
  
  
  // run the fully parallel sampler
  graphlab::timer timer;
  timer.start();
  parallel_sample(factor_graph, mrf_graph, 
                  clopts.ncpus,
                  runtime,
                  treesize,
                  priorities);
  double actual_runtime = timer.current_time();
  std::cout << "Runtime: " << actual_runtime << std::endl;

  std::cout << "Computing unnormalized log-likelihood" << std::endl;
  double loglik = unnormalized_loglikelihood(mrf_graph,
                                             factor_graph.factors());

  std::cout << "LogLikelihood: " << loglik << std::endl;
  std::cout << "Saving final prediction" << std::endl;



  //  // Plot the final answer
  // size_t rows = std::sqrt(mrf_graph.num_vertices());
  // image img(rows, rows);
  // std::vector<double> values(1);
  // for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {
  //   mrf::vertex_data& vdata = mrf_graph.vertex_data(vid);
  //   vdata.belief.normalize();
  //   vdata.belief.expectation(values);
  //   img.pixel(vid) = values[0];
  // }
  // img.save("final_pred.pgm");

  // for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {   
  //   img.pixel(vid) = mrf_graph.vertex_data(vid).updates;
  // }
  // img.save("sample_count.pgm");

  // for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {   
  //   img.pixel(vid) = mrf_graph.vertex_data(vid).updates == 0;
  // }
  // img.save("unsampled.pgm");


  

  




  return EXIT_SUCCESS;
} // end of main













#include <graphlab/macros_undef.hpp>


