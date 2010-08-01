/**
 *
 * Parallel blocked gibbs using graphlab
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <ctime>
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
#include "util.hpp"
#include "data_structures.hpp"
#include "update_functions.hpp"
#include "image_denoise.hpp"
#include "drawing.hpp"

// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>






// Command Line Parsing =======================================================>


namespace boost {
  template<>
  std::string  lexical_cast< std::string >(const std::vector<size_t>& vec) {
    std::stringstream strm;
    strm << "{" ;
    for(size_t i = 0; i < vec.size(); ++i) {
      strm << vec[i];
      if(i < vec.size() - 1) strm << ", ";
    }
    strm << "}";
    return strm.str();
  }
};




struct options {
  size_t ncpus;

  size_t num_rings;
  size_t rows;
  size_t cols;
  double sigma;
  double lambda;

  size_t tree_height;
  size_t nroots;
  bool use_weights;
  double pruning;

  std::vector<size_t> nsamples;
  std::vector<size_t> times;

  std::string drawing;
  std::string corruption;
  
  std::string problem;  
  std::string scheduler;
  std::string smoothing;
  std::string experiment;
};


std::string tree_results_fn = "tree_results.tsv";
std::string async_results_fn = "async_results.tsv";
std::string colored_results_fn = "colored_results.tsv";


size_t get_next_experiment_id(const std::string& experiment_file) {
  std::ifstream fin(experiment_file.c_str());
  size_t lines = 0;
  std::string line;
  while(getline(fin, line)) lines++;
  return lines;
}












void run_colored_samples(denoise_problem& problem,
                         const std::vector<size_t>& nsamples,
                         size_t ncpus) {
  // Setup the problem
  gl::graph& graph = problem.graph;
  factorized_model& model = problem.model;

  
  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine =
    make_colored_engine(graph, sdm, ncpus);
  assert(engine != NULL); 
  add_factors_to_sdm(sdm, model.factors());
  

  sdm.set_sync(NSAMPLES_ID,
               gl::sync_ops::sum< size_t, get_nsamples >,
               gl::apply_ops::identity_print< size_t >,
               size_t(0),
               1000);
  engine->add_terminator(nsamples_terminator);

  graphlab::timer time;
  double runtime = 0;
  // Run the experiments
  for(size_t i = 0; i < nsamples.size(); ++i) {
    size_t samples = nsamples[i] * problem.rows * problem.cols;
    sdm.set_constant(MAX_NSAMPLES_ID,  graphlab::any( samples ) );

    size_t experiment_id = get_next_experiment_id(colored_results_fn);
    std::cout << "Running Colored sampler experiment: " 
              << experiment_id << std::endl;
    std::cout << "With samples: " << samples << std::endl;

    // Run the engine
    time.start();
    engine->start();
    runtime += time.current_time();

    // Get the number of samples
    sdm.sync(graph, NSAMPLES_ID);
    size_t actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = engine->last_update_count();

    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Samples requested: " << samples << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    std::cout << "Drawing the beliefs image" << std::endl;
    draw_beliefs(graph,
                 make_filename("colored_pred_", ".pgm", experiment_id),
                 problem.rows,
                 problem.cols);
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("colored_blfs_", ".blfs", experiment_id));


    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("colored_asg_", ".asg", experiment_id));


    std::cout << "Drawing Assignment Image." << std::endl;
    draw_asg(graph,
             make_filename("colored_asg_", ".pgm", experiment_id),
             problem.rows,
             problem.cols);


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;
    
    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(graph, min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;
    
    std::ofstream fout(colored_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << problem.sigma << '\t'
         << problem.lambda << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}









void run_colored_times(denoise_problem& problem,
                       const std::vector<size_t>& times,
                       size_t ncpus) {
  // Setup the problem
  gl::graph& graph = problem.graph;
  factorized_model& model = problem.model;


  
  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine =
    make_colored_engine(graph, sdm, ncpus);
  assert(engine != NULL);
  add_factors_to_sdm(sdm, model.factors());
  
  sdm.set_sync(NSAMPLES_ID,
               gl::sync_ops::sum< size_t, get_nsamples >,
               gl::apply_ops::identity_print< size_t >,
               size_t(0),
               -1);

  graphlab::timer timer;
  double runtime = 0;
  double last_time = 0;
  // Run the experiments
  for(size_t i = 0; i < times.size(); ++i) {
    size_t timeout = times[i] - last_time;
    last_time = times[i];
    if(timeout == 0) continue;
    assert(timeout <= times[i]);
    std::cout << "Timeout: " << timeout << std::endl;
    engine->set_timeout(timeout);

    size_t experiment_id = get_next_experiment_id(colored_results_fn);
    std::cout << "Running Colored sampler Time experiment:" 
              << experiment_id << std::endl;
    std::cout << "With runtime: " << times[i] << std::endl;

    // Run the engine
    timer.start();
    engine->start();
    runtime += timer.current_time();

    // Get the number of samples
    sdm.sync(graph, NSAMPLES_ID);
    size_t actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = engine->last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    std::cout << "Drawing the beliefs image" << std::endl;
    draw_beliefs(graph,
                 make_filename("colored_pred_", ".pgm", experiment_id),
                 problem.rows,
                 problem.cols);
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("colored_blfs_", ".blfs", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("colored_asg_", ".asg", experiment_id));


    std::cout << "Drawing Assignment Image." << std::endl;
    draw_asg(graph,
             make_filename("colored_asg_", ".pgm", experiment_id),
             problem.rows,
             problem.cols);


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(graph, min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    std::ofstream fout(colored_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << problem.sigma << '\t'
         << problem.lambda << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}













void run_async_samples(denoise_problem& problem,
                       const std::vector<size_t>& nsamples,
                       size_t ncpus,
                       const std::string& scheduler) {
  // Setup the problem
  gl::graph& graph = problem.graph;
  factorized_model& model = problem.model;

  
  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine = make_engine(graph,
                                    sdm,                                   
                                    ncpus,
                                    scheduler);
  assert(engine != NULL);
  add_factors_to_sdm(sdm, model.factors());

  
  sdm.set_sync(NSAMPLES_ID,
               gl::sync_ops::sum< size_t, get_nsamples >,
               gl::apply_ops::identity_print< size_t >,
               size_t(0),
               1000);
  engine->add_terminator(nsamples_terminator);


  // Set the initial schedule
  std::vector<vertex_id_t> vertices(graph.num_vertices(), 0);
  for(vertex_id_t i = 0; i < vertices.size(); ++i) vertices[i] = i;
  std::random_shuffle(vertices.begin(), vertices.end());
  const double residual = 1.0;
  const bool UseCallback = true;
  for(size_t i = 0; i < vertices.size(); ++i) {
    gl::update_task task(vertices[i], single_sample_update<UseCallback>);
    engine->get_scheduler().add_task(task, residual);
  }

  graphlab::timer time;
  double runtime = 0;
  // Run the experiments
  for(size_t i = 0; i < nsamples.size(); ++i) {
    size_t samples = nsamples[i] * problem.rows * problem.cols;
    sdm.set_constant(MAX_NSAMPLES_ID,  graphlab::any( samples ) );

    size_t experiment_id = get_next_experiment_id(async_results_fn);
    std::cout << "Running Asynchronous sampler experiment: " 
              << experiment_id << std::endl;
    std::cout << "With samples: " << samples << std::endl;

    // Run the engine
    time.start();
    engine->start();
    runtime += time.current_time();

    // Get the number of samples
    sdm.sync(graph, NSAMPLES_ID);
    size_t actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = engine->last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Samples requested: " << samples << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    std::cout << "Drawing the beliefs image" << std::endl;
    draw_beliefs(graph,
                 make_filename("async_pred_", ".pgm", experiment_id),
                 problem.rows,
                 problem.cols);
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("async_blfs_", ".blfs", experiment_id));


    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("async_asg_", ".asg", experiment_id));


    std::cout << "Drawing Assignment Image." << std::endl;
    draw_asg(graph,
             make_filename("async_asg_", ".pgm", experiment_id),
             problem.rows,
             problem.cols);


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(graph, min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    
    std::ofstream fout(async_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << problem.sigma << '\t'
         << problem.lambda << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}













void run_async_times(denoise_problem& problem,
                     const std::vector<size_t>& times,
                     size_t ncpus,
                     const std::string& scheduler) {
  // Setup the problem
  gl::graph& graph = problem.graph;
  factorized_model& model = problem.model;

  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine = make_engine(graph,
                                    sdm,
                                    ncpus,
                                    scheduler);
  assert(engine != NULL); 
  add_factors_to_sdm(sdm, model.factors());
  

  sdm.set_sync(NSAMPLES_ID,
               gl::sync_ops::sum< size_t, get_nsamples >,
               gl::apply_ops::identity_print< size_t >,
               size_t(0),
               -1);

  // Set the initial schedule
  std::vector<vertex_id_t> vertices(graph.num_vertices(), 0);
  for(vertex_id_t i = 0; i < vertices.size(); ++i) vertices[i] = i;
  std::random_shuffle(vertices.begin(), vertices.end());
  const double residual = 1.0;
  const bool UseCallback = true;
  for(size_t i = 0; i < vertices.size(); ++i) {
    gl::update_task task(vertices[i], single_sample_update<UseCallback>);
    engine->get_scheduler().add_task(task, residual);
  }

  graphlab::timer timer;
  double runtime = 0;
  double last_time = 0;
  // Run the experiments
  for(size_t i = 0; i < times.size(); ++i) {
    size_t timeout = times[i] - last_time;
    last_time = times[i];
    if(timeout == 0) continue;
    assert(timeout <= times[i]);
    std::cout << "Timeout: " << timeout << std::endl;
    engine->set_timeout(timeout);

    size_t experiment_id = get_next_experiment_id(async_results_fn);
    std::cout << "Running Asynchronous sampler Time experiment:" 
              << experiment_id << std::endl;
    std::cout << "With runtime: " << times[i] << std::endl;

    // Run the engine
    timer.start();
    engine->start();
    runtime += timer.current_time();

    // Get the number of samples
    sdm.sync(graph, NSAMPLES_ID);
    size_t actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = engine->last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    std::cout << "Drawing the beliefs image" << std::endl;
    draw_beliefs(graph,
                 make_filename("async_pred_", ".pgm", experiment_id),
                 problem.rows,
                 problem.cols);
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("async_blfs_", ".blfs", experiment_id));

    
    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("async_asg_", ".asg", experiment_id));



    std::cout << "Drawing Assignment Image." << std::endl;
    draw_asg(graph,
             make_filename("async_asg_", ".pgm", experiment_id),
             problem.rows,
             problem.cols);


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;
    
    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(graph, min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    
    std::ofstream fout(async_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << problem.sigma << '\t'
         << problem.lambda << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}






















void run_tree_samples(denoise_problem& problem,
                      std::vector<size_t> nsamples,
                      size_t ncpus,
                      const std::string& scheduler,
                      size_t nroots,
                      size_t tree_height,
                      bool use_weights,
                      double pruning) {

  // Setup the problem
  gl::graph& graph = problem.graph;
  factorized_model& model = problem.model;
  
  // Setup the engine
  gl::thread_shared_data sdm;
  set_tree_sampler_constants(sdm,
                             nroots,
                             tree_height,
                             use_weights,
                             pruning);
  gl::iengine* engine = make_engine(graph,
                                    sdm,
                                    ncpus,
                                    scheduler);
  assert(engine != NULL);

  add_factors_to_sdm(sdm, model.factors());

  
  sdm.set_sync(NSAMPLES_ID,
               gl::sync_ops::sum< size_t, get_nsamples >,
               gl::apply_ops::identity_print< size_t >,
               size_t(0),
               1000);
  engine->add_terminator(nsamples_terminator);


  // Create the initial set of tasks
  for(size_t j = 0; j < nroots; ++j) {
    // Add the root
    vertex_id_t root = get_next_root(j, nroots, graph.num_vertices());
    // Create a task to grow a tree from the root
    gl::update_task task(root, grow_root_update);
    engine->get_scheduler().add_task(task, grow_root_residual);
  }

  graphlab::timer time;
  double runtime = 0;
  // Run the experiments
  for(size_t i = 0; i < nsamples.size(); ++i) {
    size_t samples = nsamples[i] * problem.rows * problem.cols;
    sdm.set_constant(MAX_NSAMPLES_ID,  graphlab::any( samples ) );

    size_t experiment_id = get_next_experiment_id(tree_results_fn);
    std::cout << "Running Tree sampler Samples experiment: " 
              << experiment_id << std::endl;
    std::cout << "With samples: " << samples << std::endl;

    // Run the engine
    time.start();
    engine->start();
    runtime += time.current_time();

    // Get the number of samples
    sdm.sync(graph, NSAMPLES_ID);
    size_t actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = engine->last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Samples requested: " << samples << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    std::cout << "Drawing the beliefs image" << std::endl;
    draw_beliefs(graph,
                 make_filename("tree_pred_", ".pgm", experiment_id),
                 problem.rows,
                 problem.cols);
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("tree_blfs_", ".blfs", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("tree_asg_", ".asg", experiment_id));

    std::cout << "Drawing Assignment Image." << std::endl;
    draw_asg(graph,
             make_filename("tree_asg_", ".pgm", experiment_id),
             problem.rows,
             problem.cols);



    std::cout << "Rendering Edge Weights." << std::endl;
    draw_edge_weights(graph,
                      make_filename("tree_weights_", ".pgm", experiment_id),
                      problem.rows,  problem.cols);
    
    std::cout << "Rendering Final Tree." << std::endl;
    draw_tree(graph, 
              make_filename("final_tree_", ".pgm", experiment_id),
              problem.rows, problem.cols);

    std::cout << "Rendering Final STate." << std::endl;
    draw_state(graph, 
              make_filename("final_state_", ".pgm", experiment_id),
              problem.rows, problem.cols);

    

    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(graph, min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    
    std::ofstream fout(tree_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << problem.sigma << '\t'
         << problem.lambda << '\t'
         << nroots << '\t'       
         << use_weights << '\t'
         << pruning << '\t'
         << tree_height << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}










void run_tree_times(denoise_problem& problem,
                    std::vector<size_t> times,
                    size_t ncpus,
                    const std::string& scheduler,
                    size_t nroots,
                    size_t tree_height,
                    bool use_weights,
                    double pruning) {
  
  // Setup the problem
  gl::graph& graph = problem.graph;
  factorized_model& model = problem.model;
   
  // Setup the engine
  gl::thread_shared_data sdm;
  set_tree_sampler_constants(sdm,
                             nroots,
                             tree_height,
                             use_weights,
                             pruning);
  gl::iengine* engine = make_engine(graph,
                                    sdm,
                                    ncpus,
                                    scheduler);
  assert(engine != NULL);

  add_factors_to_sdm(sdm, model.factors());
  
  sdm.set_sync(NSAMPLES_ID,
               gl::sync_ops::sum< size_t, get_nsamples >,
               gl::apply_ops::identity_print< size_t >,
               size_t(0),
               -1);

  // Create the initial set of tasks
  for(size_t j = 0; j < nroots; ++j) {
    // Add the root
    vertex_id_t root = get_next_root(j, nroots, graph.num_vertices());
    // Create a task to grow a tree from the root
    gl::update_task task(root, grow_root_update);
    engine->get_scheduler().add_task(task, grow_root_residual);
  }

  graphlab::timer timer;
  double runtime = 0;
  double last_time = 0;
  // Run the experiments
  for(size_t i = 0; i < times.size(); ++i) {
    size_t timeout = times[i] - last_time;
    last_time = times[i];
    if(timeout == 0) continue;
    assert(timeout <= times[i]);
    std::cout << "Timeout: " << timeout << std::endl;
    engine->set_timeout(timeout);

    size_t experiment_id = get_next_experiment_id(tree_results_fn);
    std::cout << "Running Tree sampler Time experiment:" 
              << experiment_id << std::endl;
    std::cout << "With runtime: " << times[i] << std::endl;

    // Run the engine
    timer.start();
    engine->start();
    runtime += timer.current_time();

    // Get the number of samples
    sdm.sync(graph, NSAMPLES_ID);
    size_t actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = engine->last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    std::cout << "Drawing the beliefs image" << std::endl;
    draw_beliefs(graph,
                 make_filename("tree_pred_", ".pgm", experiment_id),
                 problem.rows,
                 problem.cols);
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("tree_blfs_", ".blfs", experiment_id));
    
    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("tree_asg_", ".asg", experiment_id));


    std::cout << "Drawing Assignment Image." << std::endl;
    draw_asg(graph,
             make_filename("tree_asg_", ".pgm", experiment_id),
             problem.rows,
             problem.cols);


    std::cout << "Rendering Edge Weights." << std::endl;
    draw_edge_weights(graph,
                      make_filename("tree_weights_", ".pgm", experiment_id),
                      problem.rows,  problem.cols);
    
    std::cout << "Rendering Final Tree." << std::endl;
    draw_tree(graph, 
              make_filename("final_tree_", ".pgm", experiment_id),
              problem.rows, problem.cols);

    std::cout << "Rendering Final STate." << std::endl;
    draw_state(graph, 
              make_filename("final_state_", ".pgm", experiment_id),
              problem.rows, problem.cols);

    
    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(graph, min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;
    
    std::ofstream fout(tree_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << problem.sigma << '\t'
         << problem.lambda << '\t'
         << nroots << '\t'       
         << use_weights << '\t'
         << pruning << '\t'
         << tree_height << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}






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

  srand(std::time(NULL));
  
  // Parse command line arguments --------------------------------------------->
  options opts; 
  bool success = parse_command_line(argc, argv, opts);
  if(!success)  return EXIT_FAILURE;
  display_options(opts);


    denoise_problem problem;
  if(opts.problem == "") {
    // Creating image for denoising
    std::cout << "Creating a synethic image." << std::endl;
    problem = denoise_problem(opts.num_rings,
                              opts.rows,
                              opts.cols,
                              opts.sigma,
                              opts.lambda,
                              opts.smoothing,
                              opts.drawing,
                              opts.corruption);
    std::cout << "Saving image. " << std::endl;
    problem.original.save("original.pgm");
    std::cout << "Saving corrupted image. " << std::endl;
    problem.noisy.save("noisy.pgm");
  } else {
    std::cout << "Loading problem: " << opts.problem << std::endl;
    problem.load(opts.problem);
  }
  
  std::cout << problem << std::endl;


  for(size_t i = 0; i < problem.graph.num_vertices(); ++i) {
    vertex_data& vdata = problem.graph.vertex_data(i);
    vdata.asg.uniform_sample();
  }
  
  // Create synthetic images -------------------------------------------------->
  std::cout << "Running experiment" << std::endl;


  if(opts.experiment == "tree") {
    if(opts.times.empty()) {
      run_tree_samples(problem,
                       opts.nsamples ,
                       opts.ncpus,
                       opts.scheduler,
                       opts.nroots,
                       opts.tree_height,
                       opts.use_weights,
                       opts.pruning);
    } else {
      run_tree_times(problem,
                     opts.times ,
                     opts.ncpus,
                     opts.scheduler,
                     opts.nroots,
                     opts.tree_height,
                     opts.use_weights,
                     opts.pruning);
    }
  } else if(opts.experiment == "async") {
    if(opts.times.empty()) {
      run_async_samples(problem,
                        opts.nsamples ,
                        opts.ncpus,
                        opts.scheduler);
    } else {
      run_async_times(problem,
                      opts.times,
                      opts.ncpus,
                      opts.scheduler);
    }
  } else if(opts.experiment == "colored") {
    if(opts.times.empty()) {
      run_colored_samples(problem,
                          opts.nsamples ,
                          opts.ncpus);
    } else {
      run_colored_times(problem,
                        opts.times,
                        opts.ncpus);
    }
  } else {
    std::cout << "Invalid experiment: " << opts.experiment << std::endl;
    return EXIT_FAILURE;
  }  
  
  
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
    ("ncpus",  boost_po::value<size_t>(&(opts.ncpus))->default_value(2),
     "Number of cpus to use.")
    ("nsamples",  boost_po::value< std::vector<size_t> >
     (&(opts.nsamples))->default_value( std::vector<size_t>(1,100) ),
     "A vector of value [required]")
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
    
    ("treeheight", boost_po::value<size_t>(&(opts.tree_height))->
     default_value(200*200),
     "Maximum tree size for tree sampling")
    ("weights", boost_po::value<bool>(&(opts.use_weights))->
     default_value(false),
     "Use weights when constructing trees")
    ("pruning",  boost_po::value<double>(&(opts.pruning))->default_value(0),
     "Pruning threshold (== 0 implies no pruning)")

    ("nroots", boost_po::value<size_t>(&(opts.nroots))->
     default_value(3),
     "number of simultaneous tree roots")
    
    ("times",
     boost_po::value< std::vector<size_t> >(&(opts.times)),
     "A vector of times")

    ("problem",
     boost_po::value<std::string>(&(opts.problem))->default_value(""),
     "A problem description file")

    ("drawing",
     boost_po::value<std::string>(&(opts.drawing))->default_value("sunset"),
     "{sunset, checkerboard}")

    ("corruption",
     boost_po::value<std::string>(&(opts.corruption))->default_value("gaussian"),
     "{gaussian, flip}")

    
    ("smoothing",
     boost_po::value<std::string>(&(opts.smoothing))->default_value("square"),
     "Options are {square, laplace}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("multiqueue_fifo"),
     "Options are {multiqueue_fifo, fifo, clustered_priority}")

    
    ("experiment",
     boost_po::value<std::string>(&(opts.experiment))->default_value("tree"),
     "Options are {tree, async}");


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
  std::cout << "ncpus:          " << opts.ncpus << std::endl
            << "num_rings:      " << opts.num_rings << std::endl
            << "rows:           " << opts.rows << std::endl
            << "cols:           " << opts.cols << std::endl
            << "sigma:          " << opts.sigma << std::endl
            << "lambda:         " << opts.lambda << std::endl
            << "treeheight:     " << opts.tree_height << std::endl
            << "nroots:         " << opts.nroots << std::endl
            << "weights:        " << opts.use_weights << std::endl
            << "pruning:        " << opts.pruning << std::endl
            << "problem:        " << opts.problem << std::endl
            << "smoothing:      " << opts.smoothing << std::endl
            << "drawing:        " << opts.drawing << std::endl
            << "corruption:     " << opts.corruption << std::endl
            << "scheduler:      " << opts.scheduler << std::endl
            << "experiment:     " << opts.experiment << std::endl
            << "times:          "
            << boost::lexical_cast<std::string>(opts.times) << std::endl
            << "nsamples:       "
            << boost::lexical_cast<std::string>(opts.nsamples) << std::endl;
}




