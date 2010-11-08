/**
 *
 * Parallel blocked gibbs using graphlab
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

#include "data_structures.hpp"
#include "protein_side_chain.hpp"

#include "update_functions.hpp"

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
  size_t tree_height;
  size_t nroots;
  bool use_weights;
  double pruning;

  std::vector<size_t> nsamples;
  std::vector<size_t> times;

  
  std::string protein;  
  std::string scheduler;
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







void run_colored_samples(protein_data& protein,
                         const std::vector<size_t>& nsamples,
                         size_t ncpus) {

  std::cout << "Building graph." << std::endl;
  gl::graph graph;
  construct_protein_sidechain_graph(protein, graph);  
 
  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine =
    make_colored_engine(graph, sdm, protein.edge_factors, ncpus);
  
  assert(engine != NULL); 
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
    size_t samples = nsamples[i] * protein.potentials.size();
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
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("colored_blfs_", ".blfs", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("colored_asg_", ".asg", experiment_id));


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, sdm, EDGE_FACTOR_ID);
    std::cout << "Loglikelihood:    " << loglik << std::endl;
    
    
    std::ofstream fout(colored_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << protein.potentials.size() << '\t'
         << protein.edges.size() << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}







void run_colored_times(protein_data& protein,
                       const std::vector<size_t>& times,
                       size_t ncpus) {

  std::cout << "Building graph." << std::endl;
  gl::graph graph;
  construct_protein_sidechain_graph(protein, graph);  
  
  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine =
    make_colored_engine(graph, sdm, protein.edge_factors, ncpus);
  assert(engine != NULL); 
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

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("colored_blfs_", ".blfs", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("colored_asg_", ".asg", experiment_id));


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, sdm, EDGE_FACTOR_ID);
    std::cout << "Loglikelihood:    " << loglik << std::endl;


    std::ofstream fout(colored_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << protein.potentials.size() << '\t'
         << protein.edges.size() << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}









void run_async_samples(protein_data& protein,
                       const std::vector<size_t>& nsamples,
                       size_t ncpus,
                       const std::string& scheduler) {
  // Setup the problem
  std::cout << "Building graph." << std::endl;
  gl::graph graph;
  construct_protein_sidechain_graph(protein, graph);  
  
  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine = make_engine(graph,
                                    sdm,
                                    protein.edge_factors,
                                    ncpus,
                                    scheduler);
  assert(engine != NULL); 
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
    size_t samples = nsamples[i] * protein.potentials.size();
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

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("async_blfs_", ".blfs", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("async_asg_", ".asg", experiment_id));


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, sdm, EDGE_FACTOR_ID);
    std::cout << "Loglikelihood:    " << loglik << std::endl;


    std::ofstream fout(async_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << protein.potentials.size() << '\t'
         << protein.edges.size() << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}















void run_async_times(protein_data& protein,
                     const std::vector<size_t>& times,
                     size_t ncpus,
                     const std::string& scheduler) {
  std::cout << "Building graph." << std::endl;
  gl::graph graph;
  construct_protein_sidechain_graph(protein, graph);  
  
  // Setup the engine
  gl::thread_shared_data sdm;
  gl::iengine* engine = make_engine(graph,
                                    sdm,
                                    protein.edge_factors,
                                    ncpus,
                                    scheduler);
  assert(engine != NULL); 
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
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("async_blfs_", ".blfs", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("async_asg_", ".asg", experiment_id));

    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, sdm, EDGE_FACTOR_ID);
    std::cout << "Loglikelihood:    " << loglik << std::endl;
        
    std::ofstream fout(async_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << protein.potentials.size() << '\t'
         << protein.edges.size() << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}












void run_tree_samples(protein_data& protein,
                      std::vector<size_t> nsamples,
                      size_t ncpus,
                      const std::string& scheduler,
                      size_t nroots,
                      size_t tree_height,
                      bool use_weights,
                      double pruning) {
  std::cout << "Building graph." << std::endl;
  gl::graph graph;
  construct_protein_sidechain_graph(protein, graph);  
  
  // Setup the engine
  gl::thread_shared_data sdm;
  set_tree_sampler_constants(sdm,
                             nroots,
                             tree_height,
                             use_weights,
                             pruning);
  gl::iengine* engine = make_engine(graph,
                                    sdm,
                                    protein.edge_factors,
                                    ncpus,
                                    scheduler);
  assert(engine != NULL); 
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
    size_t samples = nsamples[i] * protein.potentials.size();
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

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("tree_blfs_", ".blfs", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("tree_asg_", ".asg", experiment_id));    

    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, sdm, EDGE_FACTOR_ID);
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    
    std::ofstream fout(tree_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << protein.potentials.size() << '\t'
         << protein.edges.size() << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << nroots << '\t'       
         << use_weights << '\t'
         << pruning << '\t'
         << tree_height << std::endl;
    fout.close();
  } // end of for loop 
  delete engine;
  std::cout << "Finished." << std::endl;
}

 








void run_tree_times(protein_data& protein,
                    std::vector<size_t> times,
                    size_t ncpus,
                    const std::string& scheduler,
                    size_t nroots,
                    size_t tree_height,
                    bool use_weights,
                    double pruning) {

  std::cout << "Building graph." << std::endl;
  gl::graph graph;
  construct_protein_sidechain_graph(protein, graph);  

  // Setup the engine
  gl::thread_shared_data sdm;
  set_tree_sampler_constants(sdm,
                             nroots,
                             tree_height,
                             use_weights,
                             pruning);
  gl::iengine* engine = make_engine(graph,
                                    sdm,
                                    protein.edge_factors,
                                    ncpus,
                                    scheduler);
  assert(engine != NULL); 
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

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(graph,
                 make_filename("tree_blfs_", ".blfs", experiment_id));
    
    std::cout << "Saving Assignment." << std::endl;
    save_asg(graph,
             make_filename("tree_asg_", ".asg", experiment_id));

    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(graph, sdm, EDGE_FACTOR_ID);
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    
    std::ofstream fout(tree_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << scheduler << '\t'
         << protein.potentials.size() << '\t'
         << protein.edges.size() << '\t'
         << actual_samples << '\t'
         << ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
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

  // Parse command line arguments --------------------------------------------->
  options opts; 
  bool success = parse_command_line(argc, argv, opts);
  if(!success)  return EXIT_FAILURE;
  display_options(opts);


  // Load the protein file ---------------------------------------------------->
  protein_data protein;
  std::cout << "Loading protein file." << std::endl;
  parse_protein(opts.protein, protein);
  std::cout << protein << std::endl;
  
  
  // Create synthetic images -------------------------------------------------->

  std::cout << "Running experiment" << std::endl;
  if(opts.experiment == "tree") {
    if(opts.times.empty()) {
      run_tree_samples(protein,
                       opts.nsamples ,
                       opts.ncpus,
                       opts.scheduler,
                       opts.nroots,
                       opts.tree_height,
                       opts.use_weights,
                       opts.pruning);
    } else {
      run_tree_times(protein,
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
      run_async_samples(protein,
                        opts.nsamples ,
                        opts.ncpus,
                        opts.scheduler);
    } else {
      run_async_times(protein,
                      opts.times,
                      opts.ncpus,
                      opts.scheduler);
    }
  } else if(opts.experiment == "colored") {
    if(opts.times.empty()) {
      run_colored_samples(protein,
                          opts.nsamples ,
                          opts.ncpus);
    } else {
      run_colored_times(protein,
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
    ("times",
     boost_po::value< std::vector<size_t> >(&(opts.times)),
     "A vector of times")
    ("nsamples",  boost_po::value< std::vector<size_t> >
     (&(opts.nsamples))->default_value( std::vector<size_t>(1,100) ),
     "A vector of value [required]")
    ("treeheight", boost_po::value<size_t>(&(opts.tree_height))->
     default_value(20),
     "Maximum tree size for tree sampling")
    ("weights", boost_po::value<bool>(&(opts.use_weights))->
     default_value(false),
     "Use weights when constructing trees")
    ("pruning",  boost_po::value<double>(&(opts.pruning))->default_value(0),
     "Pruning threshold (== 0 implies no pruning)")
    ("nroots", boost_po::value<size_t>(&(opts.nroots))->
     default_value(3),
     "number of simultaneous tree roots")
    ("protein",
     boost_po::value<std::string>(&(opts.protein)),
     "The protein folder")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("multiqueue_fifo"),
     "Options are {multiqueue_fifo, fifo, clustered_priority}")   
    ("experiment",
     boost_po::value<std::string>(&(opts.experiment))->default_value("tree"),
     "Options are {tree, async}");


  // Parse the arguments
  boost_po::variables_map vm;
  boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
  boost_po::positional_options_description pos_opts;
  pos_opts.add("protein",1);
  store(boost_po::command_line_parser(argc, argv)
        .options(desc).positional(pos_opts).run(), vm);
  boost_po::notify(vm);
  if(vm.count("help") || !vm.count("protein") ) {
    std::cout << desc << std::endl;
    return false;
  }
  return true;
} // end of parse command line arguments


void display_options(options& opts) {
  std::cout << "ncpus:          " << opts.ncpus << std::endl
            << "tree_height:    " << opts.tree_height << std::endl
            << "nroots:         " << opts.nroots << std::endl
            << "weights:        " << opts.use_weights << std::endl
            << "pruning:        " << opts.pruning << std::endl
            << "protein:        " << opts.protein << std::endl
            << "scheduler:      " << opts.scheduler << std::endl
            << "experiment:     " << opts.experiment << std::endl
            << "times:          "
            << boost::lexical_cast<std::string>(opts.times) << std::endl
            << "nsamples:       "
            << boost::lexical_cast<std::string>(opts.nsamples) << std::endl;

}




