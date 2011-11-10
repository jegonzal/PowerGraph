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


/*
 *  PageRank tutorial application.
 *  pagerank.cpp
 *
 *  For description of the PageRank algorithm, see Wikipedia article
 *  http://en.wikipedia.org/wiki/Pagerank
 */

#include <iomanip>

#include "pagerank.hpp"



#include <graphlab/macros_def.hpp>

//! Global random reset probability
double RANDOM_RESET_PROBABILITY = 0.15;

//! Global accuracy tolerance
double ACCURACY = 1e-5;

/**
 * The factorized page rank update function
 */
class pagerank_update : 
  public graphlab::iupdate_functor<graph_type, pagerank_update> {
  float prio;
public:
  pagerank_update(const float& prio = 0) : prio(prio) { }
  double priority() const { return prio; }
  void operator+=(const pagerank_update& other) { prio += other.prio; }
  void operator()(icontext_type& context) {
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    // Compute weighted sum of neighbors
    float sum = vdata.value * vdata.self_weight;    
    foreach(edge_id_type eid, context.in_edge_ids()) 
      sum += context.edge_data(eid).weight * 
        context.neighbor_vertex_data(context.source(eid)).value;
    // Add random reset probability
    sum = RANDOM_RESET_PROBABILITY/context.num_vertices() + 
      (1-RANDOM_RESET_PROBABILITY)*sum;
    vdata.old_value = vdata.value;
    vdata.value = sum;
    foreach(edge_id_type eid, context.out_edge_ids()) {    
      const float residual = context.edge_data(eid).weight * 
        std::fabs(vdata.value - vdata.old_value);
      // If the neighbor changed sufficiently add to scheduler.
      if(residual > ACCURACY) 
        context.schedule(context.target(eid), pagerank_update(residual));      
    }
  } // end of operator()  
}; // end of pagerank update functor




int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_file;
  std::string format = "metis";
  std::string update_type = "basic";
  clopts.attach_option("graph",
                       &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.attach_option("format",
                       &format, format,
                       "The graph file format: {metis, jure, tsv}");
  clopts.attach_option("accuracy",
                       &ACCURACY, ACCURACY,
                       "residual termination threshold");
  clopts.attach_option("resetprob",
                       &RANDOM_RESET_PROBABILITY, RANDOM_RESET_PROBABILITY,
                       "Random reset probability");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Termination bound:  " << ACCURACY 
            << std::endl
            << "Reset probability:  " << RANDOM_RESET_PROBABILITY
            << std::endl;  
  // Setup the GraphLab execution core and load graph -------------------------
  graphlab::core<graph_type, pagerank_update> core;
  core.set_options(clopts); // attach the command line options to the core
  const bool success = load_graph(graph_file, format, core.graph());
  if(!success) {
    std::cout << "Error in reading file: " << graph_file << std::endl;
  }

  // Run the PageRank ---------------------------------------------------------
  core.schedule_all(pagerank_update(1));
  const double runtime = core.start();  // Run the engine
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  std::cout << "Update Rate (updates/second): " 
            << core.last_update_count() / runtime
            << std::endl;

 

  // Output Results -----------------------------------------------------------
  // Output the top 5 pages
  std::vector<graph_type::vertex_id_type> top_pages;
  get_top_pages(core.graph(), 5, top_pages);
  for(size_t i = 0; i < top_pages.size(); ++i) {
    std::cout << std::setw(10) << top_pages[i] << ":\t" << std::setw(10) 
              << core.graph().vertex_data(top_pages[i]).value << std::setw(10) 
              << core.graph().vertex_data(top_pages[i]).nupdates << std::endl;
  }

  // Write the pagerank vector
  std::cout << "Saving pagerank vector." << std::endl;
  save_pagerank("pagerank.tsv", core.graph());
  std::cout << "Finished." << std::endl;

  return EXIT_SUCCESS;
} // End of main
