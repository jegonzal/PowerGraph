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
 * @file dispatcher.cpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#include "dispatcher.hpp"
#include <graphlab/macros_def.hpp>

using namespace graphlab;

/////////////////////////////// INSTANCE MEMBERS ///////////////////////////////
// TODO
dispatcher_update::dispatcher_update(){}

dispatcher_update::dispatcher_update(const dispatcher_update& other){}

inline void dispatcher_update::operator+=(const dispatcher_update& other){}

void dispatcher_update::operator()(icontext_type& context){
  process p = process::get_process();
  p << "oh yeah...";
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////// MAIN METHOD /////////////////////////////////
int main(int argc, char** argv) {

  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  
  // Parse command line options -----------------------------------------------
  command_line_options clopts("GraphLab Dispatcher");
  std::string graph_file = "toy.tsv", updater;
  std::string format = "tsv";
  clopts.attach_option("graph", &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", &format, format,
                       "The graph file format: {metis, snap, tsv}");
  clopts.attach_option("updater", &updater, updater,
                       "The updater to execute.");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  // Setup the GraphLab execution core and load graph -------------------------
  core<graph_type, dispatcher_update> core;
  core.set_options(clopts);      // attach the command line options to the core
  std::cout << "Loading graph from file" << std::endl;
  const bool success = graph_ops::
    load_structure(graph_file, format, core.graph());
  if(!success) {
    std::cout << "Error in reading file: " << graph_file << std::endl;
  }

  // Run the Dispatcher -------------------------------------------------------
  process::set_executable(updater);
  core.schedule_all(dispatcher_update());
  std::cout << "Running dispatcher..." << std::endl;
  const double runtime = core.start();  // Run the engine
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  std::cout << "Update Rate (updates/second): " 
            << core.last_update_count() / runtime
            << std::endl;

  // Output Results -----------------------------------------------------------
  return EXIT_SUCCESS;
  
}
////////////////////////////////////////////////////////////////////////////////

#include <graphlab/macros_undef.hpp>
