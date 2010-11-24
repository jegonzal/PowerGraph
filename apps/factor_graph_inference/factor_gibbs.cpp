
/**
 * This file runs gibbs sampling an arbitrary discrete factor graphs
 * encoded in the Alchemy factor graph format.
 *
 *
 * @author: Joseph Gonzalez
 */


#include <cassert>
#include <iostream>
#include <fstream>


// Used to parse command line input
#include <boost/program_options.hpp>

// Contains the graphlab library
#include <graphlab.hpp>


// Support code for loading alchemy factor graphs
#include "factor_graph.hpp"

// Support code for this cpp file
#include "factor_gibbs.hpp"

// Shared variables (we will eliminate these later)
size_t total_colors;
size_t nsamples;

// Command Line Support
// ============================================================================>
struct options {
  size_t ncpus;
  size_t nsamples;
  bool use_planner;
  std::string network_filename;
  std::string beliefs_filename;
  std::string stats_filename;
  std::string engine;
  std::string scope;
  std::string scheduler;
  std::string visualizer;
};
bool parse_command_line(int argc, char** argv, options& opts);
void display_options(options& opts);




/** The standard main built around the majority of the schedulers */
void standard_main(const options& opts,
                   gl_types::graph& graph,
                   gl_types::ishared_data_manager& sdm,
                   double& runtime,
                   size_t& update_count);


/** The main built around the set_scheduler */
void set_scheduler_main(const options& opts,
                        gl_types::graph& graph,
                        gl_types::ishared_data_manager& sdm,
                        double& runtime,
                        size_t& update_count);


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
  if(!success) {
    return EXIT_FAILURE;
  }
  display_options(opts);

  nsamples = opts.nsamples;
  
  // Load the alchmey network ------------------------------------------------->
  std::cout << "Loading Factor Graph" << std::endl;
  factor_graph fgraph;
  fgraph.load_alchemy(opts.network_filename);
  std::cout << "Loaded: " << fgraph.variables().size() << " variables and  "
            << fgraph.factors().size() << " factors." << std::endl;  

  // Building Clique graph  --------------------------------------------------->
  std::cout << "Building Clique Graph" << std::endl;
  gl_types::graph mrf_graph;
  make_gibbs_graph(fgraph, mrf_graph);
  std::cout << "Finished with " << mrf_graph.num_edges() << " edges."
            << std::endl;

  // Coloring graph ----------------------------------------------------------->
  std::cout << "Parallel coloring graph." << std::endl;
  total_colors = parallel_graph_color(mrf_graph, opts.ncpus);
  std::cout << "Finished with " << total_colors << " colors." << std::endl;
  assert(total_colors > 0);
    
  // Setup extra parameters in shared data ------------------------------------>
  std::cout << "Setting up shared data." << std::endl;
  gl_types::thread_shared_data sdm;
  fill_shared_data(fgraph, sdm);
  std::cout << "Finished!" << std::endl;
  std::cout << "Graph has " << mrf_graph.num_vertices() << " vertices and " 
                            << mrf_graph.num_edges() << " edges"<< std::endl;
  
  // Create the engine -------------------------------------------------------->
  double runtime = 0;
  size_t update_count = 0;
  if(opts.scheduler == "set") {
    set_scheduler_main(opts, mrf_graph, sdm, runtime, update_count);
  } else {
    standard_main(opts, mrf_graph, sdm, runtime, update_count);
  }

  // Saving beliefs file ------------------------------------------------------>
  std::cout << "Saving beliefs." << std::endl;
  save_beliefs(opts.beliefs_filename, fgraph, mrf_graph);
  std::cout << "Finished." << std::endl;


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
        << mrf_graph.num_vertices() << ", "
        << mrf_graph.num_edges() << ", "
        << fgraph.variables().size() << ", "
        << fgraph.factors().size() << ", "
        << (opts.use_planner? 1 : 0) << ", "
        << total_colors << ", "
        << opts.nsamples << std::endl;
  stats.close();
  std::cout << "Finished saving statistics file." << std::endl;

  

  return EXIT_SUCCESS;
}  // End of main




// Implementation
// ============================================================================>
void standard_main(const options& opts,
                   gl_types::graph& graph,
                   gl_types::ishared_data_manager& sdm,
                   double& runtime,
                   size_t& update_count) {
  std::cout << "Running standard main" << std::endl;
  gl_types::iengine* engine =
    graphlab::engine_factory::new_engine(opts.engine,
                                         opts.scheduler,
                                         opts.scope,
                                         graph,
                                         opts.ncpus);
  if(engine == NULL) {
    std::cout << "Unable to construct engine!" << std::endl;
    exit( EXIT_FAILURE );   
  }
  // Set the shared data manager
  engine->set_shared_data_manager(&sdm);

  const bool use_callback = true;
  gl_types::update_function update_function =
    gibbs_update<use_callback>;

  // Set the update function
  engine->sched_options().add_option("update_function", update_function);

  engine->add_task_to_all(update_function, DEFAULT_RESIDUAL);

  // Running the engine ------------------------------------------------------->
  std::cout << "Running the engine. " << std::endl;
  graphlab::timer timer; timer.start(); 

  engine->start();
  
  runtime = timer.current_time();
  update_count = engine->last_update_count();

  
  // Display the output
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;

  std::cout << "Done!" << std::endl;  
  if(engine != NULL) delete engine;

} // standard main




/**
 * This function is called when the set scheduler is selected and uses
 * a more sophisticated setup built around the graph coloring.
 */
void set_scheduler_main(const options& opts, 
                        gl_types::graph& graph,
                        gl_types::ishared_data_manager& sdm,
                        double& runtime,
                        size_t& update_count) {
  
  std::cout << "Running set scheduler main" << std::endl;
  // Here we build a set scheduler and disable locking on the scopes
  gl_types::iengine* engine =
    graphlab::engine_factory::new_engine(opts.engine,
                                         "set",
                                         "unsync",
                                         graph,
                                         opts.ncpus);

  if(engine == NULL) {
    std::cout << "Unable to construct engine!" << std::endl;
    exit( EXIT_FAILURE );   
  }

  // Set the shared data 
  engine->set_shared_data_manager(&sdm);
  // TODO set scheduler disabled
/*  if (opts.use_planner) {
    std::cout << "Using set schedule planner." << std::endl;
     engine->get_scheduler().set_option(gl_types::scheduler_options::SCHEDULING_FUNCTION, 
                                          (void*)planned_color_schedule);
  } else {
    std::cout << "Using the basic schedule planner." << std::endl;
    engine->get_scheduler().set_option(gl_types::scheduler_options::SCHEDULING_FUNCTION, 
                                        (void*)basic_color_schedule);
  }*/


  // Running the engine ------------------------------------------------------->
  std::cout << "Running the engine. " << std::endl;
  graphlab::timer timer; timer.start();
  
  engine->start();
  
  runtime = timer.current_time();
  update_count = engine->last_update_count();

  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;

  std::cout << "Done!" << std::endl;  

  if(engine != NULL) delete engine;
}









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
    ("nsamples",  boost_po::value<size_t>(&(opts.nsamples))->default_value(100),
     "Number of samples to draw.")
    ("network",
     boost_po::value<std::string>(&(opts.network_filename)),
     "Alchemy format network file.")
    ("beliefs",
     boost_po::value<std::string>(&(opts.beliefs_filename))
     ->default_value("gibbs_beliefs.csv"),
     "Output beliefs file.")
    ("stats",
     boost_po::value<std::string>(&(opts.stats_filename))->default_value("gibbs_stats.csv"),
     "Stats file containing performance results")
    ("planner",
     boost_po::value<bool>(&(opts.use_planner))->default_value(false),
     "Options are {unsync, locked}")
    ("engine",
     boost_po::value<std::string>(&(opts.engine))->default_value("threaded"),
     "Options are {threaded, sequential}")
    ("scope",
     boost_po::value<std::string>(&(opts.scope))->default_value("locked"),
     "Options are {unsync, locked}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("fifo"),
     "Options are {fifo, priority, sampling, splash}")
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
            << "nsamples:       " << opts.nsamples << std::endl
            << "network file:   " << opts.network_filename << std::endl
            << "beliefs file:   " << opts.beliefs_filename << std::endl
            << "stats file:     " << opts.stats_filename << std::endl
            << "engine:         " << opts.engine << std::endl
            << "scope:          " << opts.scope << std::endl
            << "scheduler:      " << opts.scheduler << std::endl
            << "visualizer:     " << opts.visualizer << std::endl;

} // end of display options

