/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


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
#include "util.hpp"
#include "factorized_model.hpp"
#include "mrf.hpp"
#include "junction_tree.hpp"

#include "chromatic_sampler.hpp"
#include "jt_splash_sampler.hpp"
#include "global_variables.hpp"


// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>


// Results files =============================================================>
const std::string chromatic_results_fn = "chromatic_results.tsv";
const std::string async_results_fn     = "async_results.tsv";
const std::string splash_results_fn    = "splash_results.tsv";
const std::string jtsplash_results_fn  = "jtsplash_results.tsv";




// Command Line Arguments  ====================================================>
std::string model_filename; 
std::string experiment_type = "chromatic";
std::vector<double> runtimes(1, 10);  
bool draw_images = false;

size_t treesize = 1000;
size_t treewidth = 3;
size_t treeheight = std::numeric_limits<size_t>::max();
size_t factorsize = std::numeric_limits<size_t>::max();
size_t subthreads = 1;
bool   priorities = false;








// MAIN =======================================================================>
int main(int argc, char** argv) { 

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  // std::srand ( graphlab::timer::usec_of_day() );
  // graphlab::random::seed();
  
  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts;

  clopts.attach_option("model",
                       &model_filename, model_filename,
                       "model file name");
  clopts.add_positional("model");

  clopts.attach_option("experiment", 
                       &experiment_type, experiment_type,
                       "the type of experiment to run "
                       "{chromatic, jtsplash}");
  clopts.add_positional("experiment");

  clopts.attach_option("runtimes", 
                       &runtimes, runtimes,
                       "total runtime in seconds");


  clopts.attach_option("draw_images", 
                       &draw_images, draw_images,
                       "draw pictures (assume sqrt(numvert) rows)");


  clopts.attach_option("treesize", 
                       &treesize, treesize,
                       "The maximum number of variables in a junction tree");

  clopts.attach_option("treewidth", 
                       &treewidth, treewidth,
                       "The maximum treewidth.");

  clopts.attach_option("treeheight", 
                       &treeheight, treeheight,
                       "The maximum height of the trees. ");

  clopts.attach_option("factorsize", 
                       &factorsize, factorsize,
                       "The maximum factorsize");

  clopts.attach_option("subthreads", 
                       &subthreads, subthreads,
                       "The number of threads to use inside each tree "
                       "(zero means not used)");

  clopts.attach_option("priorities",
                       &priorities, priorities,
                       "Use priorities?");

  // Set defaults for scope and scheduler
  clopts.set_scheduler_type("fifo");
  clopts.set_scope_type("edge");

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

  // Set the global factors
  //SHARED_FACTORS.set(factor_graph.factors());
  SHARED_FACTORS_PTR = &(factor_graph.factors());
  

  std::cout << "Building graphlab MRF GraphLab core." << std::endl;
  mrf_gl::core mrf_core;
  mrf_from_factorized_model(factor_graph, mrf_core.graph());
  mrf_core.set_engine_options(clopts);

  std::cout << "Computing coloring." << std::endl;
  size_t colors = mrf_core.graph().compute_coloring();
  std::cout << "Colors: " << colors << std::endl;
  
  // Create synthetic images -------------------------------------------------->
  if(experiment_type == "chromatic") {
    run_chromatic_sampler(mrf_core, 
                          chromatic_results_fn, 
                          runtimes,
                          draw_images);
  } if(experiment_type == "jtsplash") {
    splash_settings settings; 
    settings.ntrees = mrf_core.get_engine_options().get_ncpus();
    settings.max_tree_size = treesize;
    settings.max_tree_height = treeheight;
    settings.max_tree_width = treewidth;
    settings.max_tree_height = treeheight;
    settings.priorities = priorities;
    settings.subthreads = subthreads;
    
    run_jtsplash_sampler(mrf_core.graph(),
                         jtsplash_results_fn,
                         runtimes,
                         draw_images,
                         settings);

  } else {
    std::cout << "Invalid experiment type!" << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;  
} // End of main





















#include <graphlab/macros_undef.hpp>
