#include <cassert>

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <graphlab.hpp>

#include "factor_graph.hpp"
#include "factor_bp.hpp"




// Command Line Support
// ============================================================================>

struct options {
  size_t ncpus;
  double bound;
  double damping;
  std::string network_filename;
  std::string beliefs_filename;
  std::string stats_filename;
  std::string engine;
  std::string scope;
  std::string scheduler;
  size_t splash_size;
  std::string visualizer;
};
bool parse_command_line(int argc, char** argv, options& opts);
void display_options(options& opts);

// MAIN
// ============================================================================>
int main(int argc, char** argv) {
  std::cout << "This program solves the sidechain prediction task."
            << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

  // Parse command line arguments --------------------------------------------->
  options opts; 
  bool success = parse_command_line(argc, argv, opts);
  if(!success)  return EXIT_FAILURE;
  display_options(opts);
  
  
  // Load the factor graph from file ------------------------------------------>
  std::cout << "Loading Factor Graph in Alchemy Format" << std::endl;
  factor_graph fgraph;
  fgraph.load_alchemy(opts.network_filename);
  size_t num_variables = fgraph.variables().size();
  size_t num_factors = fgraph.factors().size();
  std::cout << "Finished!" << std::endl;

  
  // Build the BP graph from the factor graph---------------------------------->
  std::cout << "Building BP graph from the factor graph" << std::endl;
  gl_types::graph bpgraph;
  make_bp_graph(fgraph, bpgraph); 
  size_t num_vertices = bpgraph.num_vertices();
  assert(num_vertices == num_variables + num_factors);
  size_t num_edges = bpgraph.num_edges();
  std::cout << "Loaded: " << num_vertices << " vertices "
            << "and " << num_edges << " edges." << std::endl;
  std::cout << "Finished!" << std::endl;

  
  // Setup extra parameters in shared data ------------------------------------>
  gl_types::thread_shared_data sdm;
  sdm.set_constant(BOUND_ID, graphlab::any(opts.bound));
  sdm.set_constant(DAMPING_ID, graphlab::any(opts.damping));

    
  // Create the engine -------------------------------------------------------->
  gl_types::iengine* engine =
    graphlab::engine_factory::new_engine(opts.engine,
                                         opts.scheduler,
                                         opts.scope,
                                         bpgraph,
                                         opts.ncpus);
  if(engine == NULL) {
    std::cout << "Unable to construct engine!" << std::endl;
    return EXIT_FAILURE;   
  }

  // Extra engine configuration ----------------------------------------------->
  // Set the shared data manager
  engine->set_shared_data_manager(&sdm);
  // Set the splash size if necessary
  engine->sched_options().add_option("splash_size", opts.splash_size);
  // Set the update function
  engine->sched_options().add_option("update_function", bp_update);

  // Tell the scheduler that the bp_update function should be applied
  // to all vertices with priority:
  double initial_priority = 100.0;
  engine->add_task_to_all(bp_update, initial_priority);

  

  // Running the engine ------------------------------------------------------->
  std::cout << "Running the engine. " << std::endl;
  graphlab::timer timer; timer.start();  
  // Run the engine (this blocks until their are no tasks left).
  // Convergence is determined when there are no update tasks left.
  engine->start();

  double runtime = timer.current_time();
  std::cout << "Done!" << std::endl;
  
  // Print some fun facts------------------------------------------------------>
  size_t update_count = engine->last_update_count();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;

  
  
  // Save the beliefs --------------------------------------------------------->
  std::cout << "Saving the beliefs. " << std::endl;
  save_beliefs(opts.beliefs_filename, fgraph, bpgraph, num_variables);
  std::cout << "Finished saving beliefs." << std::endl;


  // Save statistics about the run -------------------------------------------->
  std::cout << "Saving the statistics file. " << std::endl;
  std::ofstream stats(opts.stats_filename.c_str(), std::ios::app);
  assert(stats.good());
  stats << opts.scheduler << ", "
        << opts.network_filename << ", "
        << opts.scope << ", "
        << opts.ncpus << ", "
        << runtime << ", "
        << update_count << ", "
        << num_vertices << ", "
        << num_edges << ", "
        << num_variables << ", "
        << num_factors << ", "
        << opts.splash_size << std::endl;
  stats.close();
  std::cout << "Finished saving statistics file." << std::endl;

  
  if(engine != NULL) delete engine;
  return EXIT_SUCCESS;
}



// Implementation
// ============================================================================>


bool parse_command_line(int argc, char** argv, options& opts) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Protein side chain prediction.  Usage: \n"
         "protein_bp [network name] [options]");
  // Set the program options
  desc.add_options()
    ("help",   "produce this help message")
    ("ncpus",  boost_po::value<size_t>(&(opts.ncpus))->default_value(2),
     "Number of cpus to use.")
    ("bound",  boost_po::value<double>(&(opts.bound))->default_value(1E-5),
     "Residual termination bound")
    ("damping",  boost_po::value<double>(&(opts.damping))->default_value(0.3),
     "Residual termination bound")
    ("network",
     boost_po::value<std::string>(&(opts.network_filename)),
     "Alchemy format network file.")
    ("beliefs",
     boost_po::value<std::string>(&(opts.beliefs_filename))
     ->default_value("bp_belief.csv"),
     "Output beliefs file.")
    ("stats",
     boost_po::value<std::string>(&(opts.stats_filename))
     ->default_value("bp_stats.csv"),
     "Stats file containing performance results")
    ("engine",
     boost_po::value<std::string>(&(opts.engine))->default_value("threaded"),
     "Options are {threaded, sequential}")
    ("scope",
     boost_po::value<std::string>(&(opts.scope))->default_value("locked"),
     "Options are {unsync, locked}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("splash"),
     "Options are {fifo, priority, sampling, splash}")
    ("splashsize",
     boost_po::value<size_t>(&(opts.splash_size))->default_value(50),
     "The size of each splash when the splash scheduler is used")
    ("visualizer",
     boost_po::value<std::string>(&(opts.visualizer))->default_value("false"),
     "Use visualizer server {true, false, console}");
  // Parse the arguments
  boost_po::variables_map vm;
  boost_po::positional_options_description pos_opts;
  pos_opts.add("network",1);
  store(boost_po::command_line_parser(argc, argv)
        .options(desc).positional(pos_opts).run(), vm);
  boost_po::notify(vm);
  if(vm.count("help") || vm.count("network") == 0) {
    std::cout << desc << std::endl;
    return false;
  }
  return true;
} // end of parse command line arguments



void display_options(options& opts) {
  std::cout << "ncpus:          " << opts.ncpus << std::endl
            << "bound:          " << opts.bound << std::endl
            << "damping:        " << opts.damping << std::endl
            << "network file:   " << opts.network_filename << std::endl
            << "beliefs file:   " << opts.beliefs_filename << std::endl
            << "stats file:     " << opts.stats_filename << std::endl
            << "engine:         " << opts.engine << std::endl
            << "scope:          " << opts.scope << std::endl
            << "scheduler:      " << opts.scheduler << std::endl
            << "Splash size:    " << opts.splash_size << std::endl
            << "visualizer:     " << opts.visualizer << std::endl;

} // end of display options

