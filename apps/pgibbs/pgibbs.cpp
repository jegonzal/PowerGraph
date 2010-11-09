/**
 *
 * This program runs the various gibbs samplers
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
#include "factorized_model.hpp"
#include "mrf.hpp"
#include "junction_tree.hpp"
#include "chromatic_update_function.hpp"





// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>


// Results files =============================================================>
const std::string chromatic_results_fn = "chromatic_results.tsv";
const std::string async_results_fn     = "async_results.tsv";
const std::string splash_results_fn    = "splash_results.tsv";
const std::string jtsplash_results_fn  = "jtsplash_results.tsv";

//! Get the number of lines in a file
size_t file_line_count(const std::string& experiment_file);

//! make a filename for the base sufic and number
std::string make_filename(const std::string& base,
                          const std::string& suffix,
                          const size_t number);



//! Statistics associated with a run
struct run_statistics {
  size_t nsamples;
  size_t nchanges;
  double loglik;
  size_t min_samples;
  size_t max_samples;
  run_statistics() :
    nsamples(0), nchanges(0), loglik(0.0),
    min_samples(std::numeric_limits<size_t>::max()), max_samples(0) { }
  run_statistics(const mrf_gl::core& core);
  void print() const {
    std::cout << "nsamples:        " << nsamples << std::endl
              << "nchanges:        " << nchanges << std::endl
              << "loglik:          " << loglik   << std::endl
              << "min_samples:     " << min_samples << std::endl
              << "max_samples:     " << max_samples << std::endl;
  }
};



//! Draw images for the problem
void draw_mrf(const size_t experiment_id,
              const std::string& base_name, 
              const mrf_graph_type& mrf);


//! Run the chromatic sampler for a fixed ammount of time
void run_chromatic_sampler(mrf_gl::core& core, 
                           const std::vector<double>& runtime);



// Command Line Arguments  ====================================================>
std::string model_filename; 
std::string experiment_type = "chromatic";
std::vector<double> runtimes(1, 10);  
bool draw_images = false;



// MAIN =======================================================================>
int main(int argc, char** argv) { 

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  std::srand ( graphlab::timer::usec_of_day() );
  graphlab::random::seed();
  
  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts;

  clopts.attach_option("model",
                       &model_filename, model_filename,
                       "model file name");
  clopts.add_positional("model");

  clopts.attach_option("experiment", 
                       &experiment_type, experiment_type,
                       "the type of experiment to run "
                       "{chromatic}");
  clopts.add_positional("experiment");

  clopts.attach_option("runtimes", 
                       &runtimes, runtimes,
                       "total runtime in seconds");

  clopts.attach_option("draw_images", 
                       &draw_images, draw_images,
                       "draw pictures (assume sqrt(numvert) rows)");


  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }


  std::cout << "Application Options" << std::endl;
  std::cout 
    << "model:          " << model_filename << std::endl
    << "experiment:     " << experiment_type << std::endl 
    << "runtime:        " 
    << boost::lexical_cast<std::string>(runtimes) << std::endl;
   
  std::cout << "Graphlab Options" << std::endl;
  clopts.print();


  // create model filename
  std::cout << "Load alchemy file." << std::endl;
  factorized_model factor_graph;
  factor_graph.load_alchemy(model_filename);


  std::cout << "Building graphlab MRF GraphLab core." << std::endl;
  mrf_gl::core mrf_core;
  mrf_core_from_factorized_model(factor_graph, mrf_core);
  mrf_core.set_engine_options(clopts);

  std::cout << "Computing coloring." << std::endl;
  size_t colors = mrf_core.graph().compute_coloring();
  std::cout << "Colors: " << colors << std::endl;
  
  // Create synthetic images -------------------------------------------------->
  if(experiment_type == "chromatic") {
    run_chromatic_sampler(mrf_core, runtimes);
  } else {
    std::cout << "Invalid experiment type!" << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;  
} // End of main





size_t file_line_count(const std::string& experiment_file) {
  std::ifstream fin(experiment_file.c_str());
  size_t lines = 0;
  std::string line;
  while(getline(fin, line)) lines++;
  fin.close();
  return lines;
}




std::string make_filename(const std::string& base,
                          const std::string& suffix,
                          const size_t number) {
  std::stringstream strm;
  strm << base
       << std::setw(10) << std::setfill('0')
       << number
       << suffix;
  std::cout << strm.str() << std::endl;
  return strm.str();
}




run_statistics::run_statistics(const mrf_gl::core& core) {
  run_statistics stats;
  // Compute the unnormalized log likelihood
  stats.loglik = unnormalized_loglikelihood(core);
  for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
    const mrf_vertex_data& vdata = core.graph().vertex_data(vid);
    stats.nsamples += vdata.nsamples;
    stats.nchanges += vdata.nchanges;
    stats.min_samples = std::min(stats.min_samples, vdata.nsamples);
    stats.max_samples = std::max(stats.max_samples, vdata.nsamples);
  } // end of for loop
} // end of compute run statistics





void run_chromatic_sampler(mrf_gl::core& core, 
                           const std::vector<double>& runtimes) {
  // Initialize scheduler
  core.set_scheduler_type("colored");
  core.set_scope_type("null");
  // Disable sched yield to prevent thread sleeping at the end of a
  // color
  core.engine().enable_sched_yield(false);
  // Use fixed update function
  core.scheduler().set_option(mrf_gl::scheduler_options::UPDATE_FUNCTION, 
                              (void*) single_gibbs_update);
  
  double total_runtime = 0;
  double actual_total_runtime = 0;
  foreach(const double experiment_runtime, runtimes) {
    total_runtime += experiment_runtime;
    // get the experiment id
    size_t experiment_id = file_line_count(chromatic_results_fn);
    std::cout << "Running chromatic sampler experiment " << experiment_id
              << " for " << experiment_runtime << " seconds." << std::endl;
    // set the termination time
    core.engine().set_timeout(experiment_runtime);
    // Run the engine
    graphlab::timer timer;
    timer.start();
    core.start();
    double actual_experiment_runtime = timer.current_time();
    actual_total_runtime += actual_experiment_runtime;
    /// ==================================================================
    // Compute final statistics of the mode
    run_statistics stats(core);
    stats.print();
    // Save the beliefs
    save_beliefs(core.graph(),
                 make_filename("chromatic_blfs_", ".tsv", experiment_id));
    // Save the current assignments
    save_asg(core.graph(),
             make_filename("chromatic_asg_", ".asg", experiment_id));
    // Save the experiment
    std::ofstream fout(chromatic_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << total_runtime << '\t'
         << actual_total_runtime << '\t'
         << core.engine().get_ncpus() << '\t'
         << stats.nsamples << '\t'
         << stats.nchanges << '\t'
         << stats.loglik << '\t'
         << stats.min_samples << '\t'
         << stats.max_samples << std::endl;
    fout.close();
    // Plot images if desired
    if(draw_images) draw_mrf(experiment_id, "chromatic", core.graph());
  }
} // end run_chromatic sampler





void draw_mrf(const size_t experiment_id,
              const std::string& base_name, 
              const mrf_graph_type& mrf) {
  size_t rows = std::sqrt(mrf.num_vertices());
  std::cout << "Rows: " << rows << std::endl;
  image img(rows, rows);
  std::vector<double> values(1);
  factor_t belief;
  for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {
    const mrf_vertex_data& vdata = mrf.vertex_data(vid);
    belief = vdata.belief;
    belief.normalize();
    belief.expectation(values);
    img.pixel(vid) = values[0];
  }
  img.pixel(0) = 0;
  img.pixel(1) = mrf.vertex_data(0).variable.arity-1;
  img.save(make_filename(base_name + "_pred_", ".pgm", experiment_id).c_str());
  
  for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {   
    img.pixel(vid) = mrf.vertex_data(vid).nsamples;
  }
  img.save(make_filename(base_name + "_updates_", ".pgm", experiment_id).c_str());
  
  for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {   
    img.pixel(vid) = mrf.vertex_data(vid).nsamples == 0;
  }
  img.save(make_filename(base_name + "_unsampled_", ".pgm", experiment_id).c_str());
  
  for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {   
    img.pixel(vid) = mrf.vertex_data(vid).asg;
  }
  img.pixel(0) = 0;
  img.pixel(1) = mrf.vertex_data(0).variable.arity-1;
  img.save(make_filename(base_name + "_final_sample_", ".pgm", experiment_id).c_str());
} // end of draw_mrf













#include <graphlab/macros_undef.hpp>
